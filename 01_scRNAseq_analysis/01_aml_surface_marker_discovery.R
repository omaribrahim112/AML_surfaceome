# Patient-wise discovery of AML-enriched surface genes from the scRNA-seq atlas.
# Input: AML_atlas_seurat.rds. Output: marker_discovery_results/. See README.md.

library(Seurat)
library(data.table)
library(parallel)

# ---- settings ----
atlas_file      <- "AML_atlas_seurat.rds"
out_dir         <- "marker_discovery_results"
assay           <- "SCT"
malignant_label <- "Blast/Tumor"
aml_group       <- "AML"
normal_group    <- "Normal_CD34"
exclude_labels  <- character(0)     # cell_type labels to drop, e.g. "Unknown"

# reference cells: "sample_plus_normal" (own nonmalignant + CD34 normals), "sample_only", "all_nonmalignant"
reference_scope <- "sample_plus_normal"

min_pct         <- 0.3
min_cells       <- 3        # min cells per group; also min size of a reference cell type
padj_cutoff     <- 0.05     # per-sample
frac_higher     <- 2/3
min_recurrence  <- 0.8
fdr_cutoff      <- 1e-8     # mean per-sample adjusted P
p_floor         <- 1e-300   # adjusted P of 0 is set to this before averaging; "FDR" below = mean adjusted P
n_cores         <- 4

# per-gene summary across samples; a gene missing from a comparison (below min.pct) fails it
summarise_genes <- function(de, sample_info, padj_cutoff, p_floor) {
  de <- as.data.table(de)
  n_eval <- nrow(sample_info)
  combined <- de[comparison == "combined"]
  pairwise <- de[comparison != "combined"]

  pw <- pairwise[, .(n_up_sig = sum(avg_log2FC > 0 & p_val_adj < padj_cutoff)), by = .(sample_id, gene)]
  pw <- merge(pw, sample_info[, .(sample_id, n_types)], by = "sample_id")
  pw[, higher_than_each_type := n_up_sig == n_types]

  sg <- merge(combined[, .(sample_id, gene, avg_log2FC, p_val_adj)],
              pw[, .(sample_id, gene, higher_than_each_type)], by = c("sample_id", "gene"), all = TRUE)
  sg[is.na(higher_than_each_type), higher_than_each_type := FALSE]
  sg[, higher_than_combined := !is.na(avg_log2FC) & avg_log2FC > 0]
  sg[, meets_criteria := higher_than_each_type & higher_than_combined]

  out <- sg[, .(Av_LogFC = mean(avg_log2FC, na.rm = TRUE),
                Sum_LogFC = sum(avg_log2FC, na.rm = TRUE),
                FDR = mean(pmax(p_val_adj, p_floor), na.rm = TRUE),
                n_samples_tested = sum(!is.na(avg_log2FC)),
                n_higher_than_each_type = sum(higher_than_each_type),
                n_higher_than_combined = sum(higher_than_combined),
                n_meets_criteria = sum(meets_criteria),
                samples_meeting_criteria = paste(sample_id[meets_criteria], collapse = "/")),
            by = .(Gene = gene)]
  out[, n_evaluable := n_eval]
  out[, frac_higher_than_combined := n_higher_than_combined / n_eval]
  out[, sample_freq := n_meets_criteria / n_eval]
  out
}

# ---- load atlas ----
dir.create(file.path(out_dir, "per_sample"), recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(atlas_file))
atlas <- readRDS(atlas_file)
md <- atlas@meta.data
stopifnot(c("sample_id", "sample_group", "cell_type") %in% names(md),
          malignant_label %in% md$cell_type,
          reference_scope %in% c("sample_plus_normal", "sample_only", "all_nonmalignant"))

aml_samples <- sort(unique(md$sample_id[md$sample_group %in% aml_group]))
message(length(aml_samples), " AML samples, ",
        length(unique(md$sample_id[md$sample_group %in% normal_group])), " normal samples")

# ---- DE for one AML sample (results cached in per_sample/) ----
test_sample <- function(s) {
  # settings in the cache name so old results are not reused
  tag <- paste(reference_scope, min_pct, min_cells, sep = "_")
  f_cache <- file.path(out_dir, "per_sample", paste0(s, "_", tag, ".rds"))
  if (file.exists(f_cache)) return(readRDS(f_cache))

  is_mal <- md$cell_type %in% malignant_label
  usable <- !is.na(md$cell_type) & !md$cell_type %in% exclude_labels
  in_scope <- switch(reference_scope,
                     sample_only        = md$sample_id %in% s,
                     sample_plus_normal = md$sample_id %in% s | md$sample_group %in% normal_group,
                     all_nonmalignant   = rep(TRUE, nrow(md)))
  case <- rownames(md)[md$sample_id %in% s & is_mal]
  ref <- rownames(md)[usable & !is_mal & in_scope]

  make_summary <- function(ref_types = character(0), note = "") {
    data.table(sample_id = s, n_malignant = length(case), n_reference = length(ref),
               n_types = length(ref_types), types_compared = paste(ref_types, collapse = ";"),
               evaluable = length(ref_types) > 0, note = note)
  }
  res <- list(summary = make_summary(note = "too few malignant or reference cells"), de = NULL)

  if (length(case) >= min_cells && length(ref) >= min_cells) {
    obj <- subset(atlas, cells = c(case, ref))
    DefaultAssay(obj) <- assay
    Idents(obj) <- "cell_type"
    # needed when the subset holds several SCT models
    if (length(methods::slot(obj[[assay]], "SCTModel.list")) > 1) obj <- PrepSCTFindMarkers(obj)

    n_by_type <- table(Idents(obj))
    ref_types <- setdiff(names(n_by_type)[n_by_type >= min_cells], malignant_label)

    if (length(ref_types) > 0) {
      fm <- function(ident_2, label) {
        r <- FindMarkers(obj, ident.1 = malignant_label, ident.2 = ident_2, test.use = "wilcox",
                         min.pct = min_pct, logfc.threshold = 0, min.cells.group = min_cells,
                         verbose = FALSE)
        if (nrow(r) == 0) return(NULL)
        data.table(sample_id = s, comparison = label, gene = rownames(r), r)
      }
      res$de <- rbindlist(c(list(fm(NULL, "combined")), lapply(ref_types, function(t) fm(t, t))))
      res$summary <- make_summary(ref_types)
    } else {
      res$summary <- make_summary(note = "no reference cell type with enough cells")
    }
  }
  saveRDS(res, f_cache)
  res
}

results <- mclapply(aml_samples, test_sample, mc.cores = n_cores)
if (any(sapply(results, inherits, "try-error"))) stop("DE failed for some samples, check per_sample/")

sample_summary <- rbindlist(lapply(results, `[[`, "summary"))
fwrite(sample_summary, file.path(out_dir, "sample_summary.tsv"), sep = "\t")
evaluable <- sample_summary[evaluable == TRUE]
message(nrow(evaluable), " of ", nrow(sample_summary), " AML samples evaluable")

# ---- cross-sample criteria ----
de <- rbindlist(lapply(results, `[[`, "de"))
genes <- summarise_genes(de, evaluable, padj_cutoff, p_floor)
genes[, higher_in_AML := n_higher_than_each_type >= 1 & frac_higher_than_combined >= frac_higher]
genes[, passes_fdr := !is.na(FDR) & FDR < fdr_cutoff]
genes[, passes_recurrence := sample_freq >= min_recurrence]

# ---- plasma membrane genes (GO:0005886 incl. child terms) ----
pm_genes <- unique(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = "GO:0005886",
                                         keytype = "GOALL", columns = "SYMBOL")$SYMBOL)
genes[, GO_plasma_membrane := Gene %in% pm_genes]
genes[, candidate := GO_plasma_membrane & higher_in_AML & passes_fdr & passes_recurrence]

# ---- write ----
first_cols <- c("Gene", "Av_LogFC", "FDR", "sample_freq")
setcolorder(genes, c(first_cols, setdiff(names(genes), first_cols)))
setorder(genes, -Av_LogFC)
fwrite(genes, file.path(out_dir, "gene_statistics.tsv"), sep = "\t")
fwrite(genes[candidate == TRUE], file.path(out_dir, "surface_candidates.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

surf <- genes[GO_plasma_membrane == TRUE]
message("plasma membrane genes higher in AML:  ", sum(surf$higher_in_AML))
message("  and mean adjusted P < ", fdr_cutoff, ":     ", sum(surf$higher_in_AML & surf$passes_fdr))
message("  and sample_freq >= ", min_recurrence, ":      ", sum(surf$candidate))
