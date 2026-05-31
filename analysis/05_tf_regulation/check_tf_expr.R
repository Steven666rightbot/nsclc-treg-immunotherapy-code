mat <- readRDS('D:/Research/potato/data/roosted/treg_207422.rds')
meta <- read.csv('D:/Research/tomato/data/gse207422_cellchat_input_corrected/weak_metadata.csv', stringsAsFactors=FALSE)
weak_bc <- meta[[1]]
mat_w <- mat[, colnames(mat) %in% weak_bc]

tfs <- c('MEF2C', 'TEAD1', 'NFKB1', 'STAT2', 'IRF8')
for (tf in tfs) {
  if (tf %in% rownames(mat_w)) {
    expr <- mat_w[tf, ]
    n_pos <- sum(expr > 0)
    cat(sprintf('%s: n_cells_expr=%d, mean=%.4f, max=%.1f\n', tf, n_pos, mean(expr), max(expr)))
  } else {
    cat(sprintf('%s: NOT IN DATA\n', tf))
  }
}
