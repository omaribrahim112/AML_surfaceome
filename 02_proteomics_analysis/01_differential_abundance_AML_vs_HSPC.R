# limma differential abundance, primary AML vs CD34+ HSPCs (TMT global proteome).
# Input: proteomics_global_protein_table.tsv, sc_candidate_genes.txt, surface_genes.txt
#        (and optionally proteomics_samples.tsv).
# Output: proteomics_results/. Rules and thresholds are described in README.md.

library(limma)

# ---- settings ----
protein_file  <- "proteomics_global_protein_table.tsv"
sample_file   <- "proteomics_samples.tsv"
sc_gene_file  <- "sc_candidate_genes.txt"
surface_file  <- "surface_genes.txt"   # plasma-membrane genes with a defined extracellular domain
out_dir       <- "proteomics_results"
n_expected    <- c(AML = 15, HSPC = 12)
min_per_group <- 3
adjp_max      <- 0.10
log2fc_min    <- 0.5
det_aml_min   <- 0.5
det_hspc_max  <- NA    # HSPC detection limit not applied (NA drops this criterion)

# ---- read table and pick samples ----
clean_gene <- function(x) {
  x <- sub("(_HUMAN|\\s*\\(HUMAN\\))$", "", as.character(x), ignore.case = TRUE)
  toupper(trimws(sub("^HUMAN[_.]", "", x)))
}

stopifnot(file.exists(protein_file), file.exists(sc_gene_file), file.exists(surface_file))
tab <- read.delim(protein_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
stopifnot("prot_acc" %in% names(tab))
gene_col <- intersect(c("gene", "Gene"), names(tab))
tab$Gene <- clean_gene(if (length(gene_col)) tab[[gene_col[1]]] else tab$prot_acc)
tab$score <- if ("prot_score" %in% names(tab)) suppressWarnings(as.numeric(tab$prot_score)) else NA_real_

quant <- grep("^N_log2_R[0-9A-Z]+ \\(", names(tab), value = TRUE)
sample_id <- sub("\\)$", "", sub("^N_log2_R[0-9A-Z]+ \\(", "", quant))
group <- ifelse(grepl("^AML\\.", sample_id), "AML", ifelse(grepl("^HSC\\.", sample_id), "HSPC", NA))
keep <- !is.na(group)
if (file.exists(sample_file)) keep <- keep & sample_id %in% read.delim(sample_file)$sample_id
quant <- quant[keep]; sample_id <- sample_id[keep]
group <- factor(group[keep], levels = c("HSPC", "AML"))

n_found <- table(group)
if (!all(n_found[names(n_expected)] == n_expected)) {
  stop("found ", n_found[["AML"]], " AML and ", n_found[["HSPC"]], " HSPC samples, expected ",
       n_expected[["AML"]], " and ", n_expected[["HSPC"]], "; list the samples in ", sample_file)
}

expr <- as.matrix(tab[, quant])
storage.mode(expr) <- "double"
expr[!is.finite(expr)] <- NA
rownames(expr) <- make.unique(tab$prot_acc)
colnames(expr) <- sample_id
is_aml <- group == "AML"
is_hspc <- group == "HSPC"

# ---- limma ----
n_aml <- rowSums(!is.na(expr[, is_aml]))
n_hspc <- rowSums(!is.na(expr[, is_hspc]))
tested <- n_aml >= min_per_group & n_hspc >= min_per_group
tested[tested] <- apply(expr[tested, ], 1, var, na.rm = TRUE) > 0
message(sum(tested), " of ", nrow(expr), " proteins tested")

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr[tested, ], design)
fit <- contrasts.fit(fit, makeContrasts(AML - HSPC, levels = design))
fit <- eBayes(fit)                    # default moderation, no robust option
tt <- topTable(fit, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")

res <- data.frame(prot_acc = rownames(tt), Gene = tab$Gene[match(rownames(tt), rownames(expr))],
                  score = tab$score[match(rownames(tt), rownames(expr))],
                  log2FC = tt$logFC, p = tt$P.Value, adjP = tt$adj.P.Val,
                  n_AML = n_aml[rownames(tt)], n_HSPC = n_hspc[rownames(tt)],
                  stringsAsFactors = FALSE)
res$detection_AML <- res$n_AML / n_expected[["AML"]]
res$detection_HSPC <- res$n_HSPC / n_expected[["HSPC"]]

# one protein per gene: highest protein score, then smallest adjusted P
res <- res[order(res$Gene, -res$score, res$adjP, na.last = TRUE), ]
by_gene <- res[!is.na(res$Gene) & res$Gene != "" & !duplicated(res$Gene), ]

# ---- candidate calls ----
sc_genes <- clean_gene(readLines(sc_gene_file))
sc_genes <- unique(sc_genes[nzchar(sc_genes)])
surface_genes <- clean_gene(readLines(surface_file))
surface_genes <- unique(surface_genes[nzchar(surface_genes)])

stats_ok <- by_gene$adjP <= adjp_max & by_gene$log2FC >= log2fc_min &
            by_gene$detection_AML >= det_aml_min
stats_ok[is.na(stats_ok)] <- FALSE
hspc_ok <- if (is.na(det_hspc_max)) rep(TRUE, nrow(by_gene)) else by_gene$detection_HSPC <= det_hspc_max
surface_ok <- by_gene$Gene %in% surface_genes
by_gene$class <- ifelse(by_gene$Gene %in% sc_genes, "single_cell",
                 ifelse(stats_ok & hspc_ok & surface_ok, "proteomics", "other"))
cand <- by_gene[by_gene$class != "other", ]

message(sum(by_gene$class != "single_cell" & stats_ok), " proteins pass adjP, log2FC and AML detection; ",
        sum(by_gene$class == "proteomics"), " also pass the surface filter")
missing_sc <- setdiff(sc_genes, by_gene$Gene)
if (length(missing_sc)) message("not tested: ", paste(missing_sc, collapse = ", "))
cand <- cand[order(-cand$log2FC), ]

message(sum(sc_genes %in% by_gene$Gene), " of ", length(sc_genes), " single-cell candidates tested; ",
        sum(cand$class == "single_cell" & cand$log2FC > 0), " higher in AML, ",
        sum(cand$class == "single_cell" & cand$log2FC > 0 & cand$adjP < 0.05), " with adjusted P < 0.05")
message(sum(cand$class == "proteomics"), " additional proteins meet the AML-enriched criteria")

# ---- abundance matrix of the candidates (input to script 02) ----
# tested proteins are represented as above; candidates that were not tested use their
# highest-scoring protein in the table
cand_genes <- unique(c(sc_genes[sc_genes %in% tab$Gene], cand$Gene[cand$class == "proteomics"]))
rep_row <- sapply(cand_genes, function(g) {
  if (g %in% by_gene$Gene) return(by_gene$prot_acc[by_gene$Gene == g])
  i <- which(tab$Gene == g)
  rownames(expr)[i[order(-tab$score[i], -rowSums(!is.na(expr[i, , drop = FALSE])), na.last = TRUE)][1]]
})
abund <- expr[rep_row, order(group, sample_id), drop = FALSE]
rownames(abund) <- cand_genes

# ---- write ----
dir.create(out_dir, showWarnings = FALSE)
res_out <- res[, c("Gene", "prot_acc", "score", "log2FC", "p", "adjP", "n_AML", "n_HSPC", "detection_AML", "detection_HSPC")]
write.table(res_out, file.path(out_dir, "all_tested_proteins_AML_vs_HSPC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cand[, c("Gene", "log2FC", "p", "adjP", "class", "prot_acc", "detection_AML", "detection_HSPC")],
            file.path(out_dir, "candidates_AML_vs_HSPC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Gene = rownames(abund), abund, check.names = FALSE),
            file.path(out_dir, "candidate_abundance_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
