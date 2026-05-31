library(CellChat)
for(i in 7:26) {
  f <- sprintf('D:/Research/tomato/results/cellchat_ko_corrected/bootstrap/cellchat_weak_bootstrap_%02d.rds', i)
  if(!file.exists(f)) {
    cat(sprintf('Iter %02d: FILE NOT FOUND\n', i))
    next
  }
  cc <- readRDS(f)
  df <- subsetCommunication(cc)
  ft <- df[df$source == 'Fibroblast' & df$target == 'Treg', ]
  prob <- if(nrow(ft) > 0) sum(ft$prob) else 0
  cat(sprintf('Iter %02d: n=%d prob=%.6e\n', i, nrow(ft), prob))
  rm(cc, df, ft)
  gc()
}
