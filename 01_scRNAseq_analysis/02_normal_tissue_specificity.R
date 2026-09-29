# GTEx (RNA) and HPA (protein) tissue specificity of the surface candidates from script 01.
# Input: marker_discovery_results/surface_candidates.tsv, GTEx and HPA downloads. Output: tissue_specificity_results/.
# Tests: one-sided Welch t test, BH; details in README.md.

library(data.table)
library(readr)

# ---- settings ----
candidate_file   <- "marker_discovery_results/surface_candidates.tsv"
gtex_tpm_file    <- "GTEx_gene_tpm.gct.gz"
gtex_attr_file   <- "GTEx_sample_attributes.txt"
hpa_file         <- "HPA_normal_tissue.tsv"
out_dir          <- "tissue_specificity_results"
fdr_cutoff       <- 0.05
hpa_bm_lymphoid  <- c("bone marrow", "lymph node", "spleen", "tonsil", "appendix", "thymus")
hpa_pooled_label <- "bone marrow and lymphoid"
gtex_bm_lymphoid <- c("Blood", "Spleen", "Lymphocytes")   # GTEx has no bone marrow category

# ---- functions ----

# Peirce's criterion (Ross 2003): N observations, k doubtful, m unknowns (the mean).
# Returns the squared deviation threshold in units of variance.
peirce_x2 <- function(N, k, m = 1) {
  if (N <= 1) return(0)
  Q <- (k^(k / N) * (N - k)^((N - k) / N)) / N
  r_new <- 1; r_old <- 0; x2 <- 0; iter <- 0
  while (abs(r_new - r_old) > N * 2e-16 && iter < 1000) {
    iter <- iter + 1
    ldiv <- r_new^k
    if (ldiv == 0) ldiv <- 1e-6
    lambda <- ((Q^N) / ldiv)^(1 / (N - k))
    x2 <- 1 + (N - m - k) / k * (1 - lambda^2)
    if (x2 < 0) {
      x2 <- 0; r_old <- r_new
    } else {
      r_old <- r_new
      r_new <- exp((x2 - 1) / 2) * 2 * pnorm(-sqrt(x2))
    }
  }
  x2
}

# TRUE for outliers (either tail). Doubtful observations are added one at a time
# as long as at least that many outliers are found.
peirce_outliers <- function(x) {
  N <- length(x)
  out <- rep(FALSE, N)
  s <- sd(x)
  if (N < 3 || !is.finite(s) || s == 0) return(out)
  k <- 1
  repeat {
    cand <- abs(x - mean(x)) > sqrt(peirce_x2(N, k)) * s
    if (sum(cand) < k) break
    out <- cand
    if (k >= N - 2) break
    k <- k + 1
  }
  out
}

# one-sided Welch t test, returns c(t, p); NA if a group is too small or constant
welch_greater <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(c(NA_real_, NA_real_))
  r <- tryCatch(t.test(x, y, alternative = "greater"), error = function(e) NULL)
  if (is.null(r)) return(c(NA_real_, NA_real_))
  c(unname(r$statistic), r$p.value)
}

# GTEx tissue (SMTSD) -> major category. Multi-region tissues take the text before " - "
# (Brain, Artery, Skin ...); the others are mapped below. Cell lines are dropped (NA).
gtex_major_tissue <- function(smtsd) {
  special <- c("Adrenal Gland" = "Adrenal", "Cervix" = "Vagina", "Fallopian Tube" = "Ovary",
               "Minor Salivary Gland" = "Salivary", "Small Intestine" = "Intestine",
               "Whole Blood" = "Blood",
               "Cells - EBV-transformed lymphocytes" = "Lymphocytes",
               "Cells - Transformed fibroblasts" = "Fibroblast",
               "Cells - Cultured fibroblasts" = "Fibroblast")
  base <- sub(" - .*$", "", smtsd)
  out <- ifelse(smtsd %in% names(special), special[smtsd],
                ifelse(base %in% names(special), special[base], base))
  out[startsWith(smtsd, "Cells - ") & !smtsd %in% names(special)] <- NA_character_
  unname(out)
}

# tpm: genes x samples matrix; tissue: major category of each column
gtex_specificity <- function(tpm, tissue) {
  idx <- split(seq_along(tissue), tissue)
  out <- vector("list", nrow(tpm))
  for (i in seq_len(nrow(tpm))) {
    v <- tpm[i, ]
    tissue_mean <- sapply(idx, function(j) mean(v[j]))
    high <- names(tissue_mean)[peirce_outliers(tissue_mean) & tissue_mean > mean(tissue_mean)]
    if (length(high) == 0) next
    out[[i]] <- rbindlist(lapply(high, function(t) {
      j <- idx[[t]]
      if (tissue_mean[[t]] < mean(v[-j])) return(NULL)
      w <- welch_greater(v[j], v[-j])
      data.table(gene = rownames(tpm)[i], tissue = t, n_samples = length(j),
                 mean_tpm = tissue_mean[[t]], mean_tpm_rest = mean(v[-j]),
                 fold_difference = tissue_mean[[t]] / mean(v[-j]), t_statistic = w[1], p_value = w[2])
    }))
  }
  res <- rbindlist(out)
  if (nrow(res) == 0) return(data.table(gene = character(), tissue = character(), p_value = numeric(), fdr = numeric()))
  res[, fdr := p.adjust(p_value, method = "BH")]
  res[]
}

# d: gene, tissue, level (0-3), one row per tissue / cell type entry.
# Tissues in pooled_tissues become one category; every other tissue is its own category.
hpa_specificity <- function(d, pooled_tissues, pooled_label) {
  d <- copy(as.data.table(d))
  d[, category := ifelse(tissue %in% pooled_tissues, pooled_label, tissue)]
  out <- list()
  for (g in unique(d$gene)) {
    x <- d[gene == g]
    for (cat in unique(x$category)) {
      inside <- x$level[x$category == cat]
      rest <- x$level[x$category != cat]
      w <- welch_greater(inside, rest)
      out[[length(out) + 1]] <- data.table(gene = g, tissue = cat, n_entries = length(inside),
                                           mean_level = mean(inside), mean_level_rest = mean(rest),
                                           t_statistic = w[1], p_value = w[2])
    }
  }
  res <- rbindlist(out)
  if (nrow(res) == 0) return(data.table(gene = character(), tissue = character(), p_value = numeric(), fdr = numeric()))
  res[, fdr := p.adjust(p_value, method = "BH")]
  res[]
}

# Reads only the candidate genes from a GTEx .gct(.gz) (2 header lines, then Name, Description,
# samples). For duplicated gene symbols the Ensembl entry with the highest mean TPM is kept.
read_gtex_tpm <- function(path, genes) {
  keep <- function(x, pos) x[x$Description %in% genes, ]
  d <- read_tsv_chunked(path, DataFrameCallback$new(keep), chunk_size = 2000, skip = 2, progress = FALSE,
                        col_types = cols(.default = col_double(),
                                         Name = col_character(), Description = col_character()))
  d <- as.data.table(d)
  if (nrow(d) == 0) stop("none of the candidate genes found in ", path)
  m <- as.matrix(d[, setdiff(names(d), c("Name", "Description")), with = FALSE])
  rownames(m) <- d$Description
  m <- m[order(-rowMeans(m)), , drop = FALSE]
  m[!duplicated(rownames(m)), , drop = FALSE]
}

# significant tests (fdr <= fdr_cutoff), one string per gene
significant_tissues <- function(res) {
  res[!is.na(fdr) & fdr <= fdr_cutoff, .(tissues = paste(tissue, collapse = ";")), by = gene]
}

# ---- run ----
stopifnot(all(file.exists(candidate_file, gtex_tpm_file, gtex_attr_file, hpa_file)))
dir.create(out_dir, showWarnings = FALSE)
cand <- fread(candidate_file)
genes <- cand$Gene

# GTEx
attrs <- fread(gtex_attr_file, select = c("SAMPID", "SMTSD"))
tpm <- read_gtex_tpm(gtex_tpm_file, genes)
tissue <- gtex_major_tissue(attrs$SMTSD[match(colnames(tpm), attrs$SAMPID)])
tpm <- tpm[, !is.na(tissue), drop = FALSE]
tissue <- tissue[!is.na(tissue)]
message("GTEx: ", ncol(tpm), " samples, ", length(unique(tissue)), " tissue categories")
gtex <- gtex_specificity(tpm, tissue)
fwrite(gtex, file.path(out_dir, "gtex_tissue_tests.tsv"), sep = "\t")

# HPA
hpa <- fread(hpa_file)
setnames(hpa, gsub(" ", "_", names(hpa)))
hpa <- hpa[Gene_name %in% genes]
hpa[, level := unname(c("Not detected" = 0, "Low" = 1, "Medium" = 2, "High" = 3)[Level])]
hpa <- hpa[!is.na(level)]                   # drops entries without a level, e.g. "Not representative"
hpa[, tissue := tolower(trimws(Tissue))]
hpa_res <- hpa_specificity(hpa[, .(gene = Gene_name, tissue, level)], hpa_bm_lymphoid, hpa_pooled_label)
fwrite(hpa_res, file.path(out_dir, "hpa_tissue_tests.tsv"), sep = "\t")

# review table
review <- copy(cand)
review[, GTEx_in_dataset := Gene %in% rownames(tpm)]
review[, HPA_in_dataset := Gene %in% hpa$Gene_name]
gtex_sig <- significant_tissues(gtex)
review[, GTEx_significant_tissues := gtex_sig$tissues[match(Gene, gtex_sig$gene)]]
review[, GTEx_BM_lymphoid_significant := Gene %in% gtex[!is.na(fdr) & fdr <= fdr_cutoff &
                                                        tissue %in% gtex_bm_lymphoid, gene]]
hpa_sig <- significant_tissues(hpa_res)
review[, HPA_significant_tissues := hpa_sig$tissues[match(Gene, hpa_sig$gene)]]
hpa_bm <- hpa_res[tissue == hpa_pooled_label]
review[, HPA_BM_lymphoid_FDR := hpa_bm$fdr[match(Gene, hpa_bm$gene)]]
review[, HPA_BM_lymphoid_significant := !is.na(HPA_BM_lymphoid_FDR) & HPA_BM_lymphoid_FDR <= fdr_cutoff]
setorder(review, -Av_LogFC)
fwrite(review, file.path(out_dir, "tissue_specificity_review.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

message(nrow(review), " candidates; ", sum(review$HPA_BM_lymphoid_significant),
        " significant in HPA bone marrow / lymphoid")
