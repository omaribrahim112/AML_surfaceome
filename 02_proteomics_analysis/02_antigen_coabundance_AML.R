# Pairwise Spearman co-abundance of the candidate antigens across the AML samples.
# Input: proteomics_results/candidate_abundance_matrix.tsv. Output: proteomics_results/antigen_coabundance_AML.tsv.

# ---- settings ----
matrix_file <- "proteomics_results/candidate_abundance_matrix.tsv"
out_file    <- "proteomics_results/antigen_coabundance_AML.tsv"
min_pairs   <- 3     # minimum samples with both proteins quantified
n_aml       <- 15

stopifnot(file.exists(matrix_file))
m <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
aml <- as.matrix(m[, grepl("^AML\\.", names(m))])
rownames(aml) <- m$Gene
if (ncol(aml) != n_aml) stop("found ", ncol(aml), " AML samples, expected ", n_aml, call. = FALSE)
message(nrow(aml), " candidates, ", ncol(aml), " AML samples")

genes <- sort(rownames(aml))
pairs <- t(combn(genes, 2))
res <- data.frame(protein1 = pairs[, 1], protein2 = pairs[, 2],
                  spearman_rho = NA_real_, p_value = NA_real_, n_pairs = NA_integer_)

for (i in seq_len(nrow(res))) {
  x <- aml[res$protein1[i], ]
  y <- aml[res$protein2[i], ]
  ok <- !is.na(x) & !is.na(y)
  res$n_pairs[i] <- sum(ok)
  if (sum(ok) >= min_pairs) {
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))   # exact P unless ties
    res$spearman_rho[i] <- unname(ct$estimate)
    res$p_value[i] <- ct$p.value
  }
}
res$p_adj_fdr <- p.adjust(res$p_value, method = "BH")
res <- res[, c("protein1", "protein2", "spearman_rho", "p_value", "p_adj_fdr", "n_pairs")]

write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(nrow(res), " pairs; ", sum(res$p_adj_fdr < 0.05, na.rm = TRUE), " with FDR < 0.05 (",
        sum(res$p_adj_fdr < 0.05 & res$spearman_rho > 0, na.rm = TRUE), " positive)")
