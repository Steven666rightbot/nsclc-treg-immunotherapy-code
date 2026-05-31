library(CellChat)

OUT <- 'D:/Research/tomato/results/cellchat_ko_corrected'

# Load weak_KO_ALL_ECM RDS and write its CSV
cat("Loading weak_KO_ALL_ECM RDS...\n")
cc <- readRDS(file.path(OUT, 'cellchat_weak_KO_ALL_ECM.rds'))
df <- subsetCommunication(cc)
write.csv(df, file.path(OUT, 'interactions_weak_KO_ALL_ECM.csv'), row.names=FALSE)
cat('Wrote interactions_weak_KO_ALL_ECM.csv\n')

# Build full summary from all existing objects
ko_list <- list(
  KO_COL1A1 = c('COL1A1'),
  KO_COL1A2 = c('COL1A2'),
  KO_COL1A1_COL1A2 = c('COL1A1','COL1A2'),
  KO_FN1 = c('FN1'),
  KO_ALL_ECM = c('COL1A1','COL1A2','FN1')
)

# Load base objects for reference
cat("Loading base objects...\n")
base_strong <- readRDS(file.path(OUT, 'cellchat_strong.rds'))
base_weak <- readRDS(file.path(OUT, 'cellchat_weak.rds'))
ft_strong_base <- subsetCommunication(base_strong)
ft_strong_base <- ft_strong_base[ft_strong_base$source=='Fibroblast' & ft_strong_base$target=='Treg',]
ft_weak_base <- subsetCommunication(base_weak)
ft_weak_base <- ft_weak_base[ft_weak_base$source=='Fibroblast' & ft_weak_base$target=='Treg',]

base_s <- sum(ft_strong_base$prob)
base_w <- sum(ft_weak_base$prob)
cat(sprintf("Base Strong Fib->Treg prob: %.6e (%d interactions)\n", base_s, nrow(ft_strong_base)))
cat(sprintf("Base Weak Fib->Treg prob: %.6e (%d interactions)\n", base_w, nrow(ft_weak_base)))

summary <- data.frame()

for (ko_name in names(ko_list)) {
  # Strong
  rds_s <- file.path(OUT, paste0('cellchat_strong_', ko_name, '.rds'))
  if (file.exists(rds_s)) {
    cat(sprintf("[%s Strong] loading...\n", ko_name))
    cc_s <- readRDS(rds_s)
    ft_s <- subsetCommunication(cc_s)
    ft_s <- ft_s[ft_s$source=='Fibroblast' & ft_s$target=='Treg',]
    ko_s <- if(nrow(ft_s)>0) sum(ft_s$prob) else 0
    retained_s <- if(base_s>0) round(ko_s/base_s*100, 2) else NA
    lost_s <- if(base_s>0) round(100 - retained_s, 2) else NA
    summary <- rbind(summary, data.frame(
      KO=ko_name, Group='Strong',
      Base_Interactions=nrow(ft_strong_base), KO_Interactions=nrow(ft_s),
      Base_Prob=base_s, KO_Prob=ko_s,
      Retained_Pct=retained_s, Lost_Pct=lost_s
    ))
  }
  # Weak
  rds_w <- file.path(OUT, paste0('cellchat_weak_', ko_name, '.rds'))
  if (file.exists(rds_w)) {
    cat(sprintf("[%s Weak] loading...\n", ko_name))
    cc_w <- readRDS(rds_w)
    ft_w <- subsetCommunication(cc_w)
    ft_w <- ft_w[ft_w$source=='Fibroblast' & ft_w$target=='Treg',]
    ko_w <- if(nrow(ft_w)>0) sum(ft_w$prob) else 0
    retained_w <- if(base_w>0) round(ko_w/base_w*100, 2) else NA
    lost_w <- if(base_w>0) round(100 - retained_w, 2) else NA
    summary <- rbind(summary, data.frame(
      KO=ko_name, Group='Weak',
      Base_Interactions=nrow(ft_weak_base), KO_Interactions=nrow(ft_w),
      Base_Prob=base_w, KO_Prob=ko_w,
      Retained_Pct=retained_w, Lost_Pct=lost_w
    ))
  }
}

cat("\n========== KO SUMMARY ==========\n")
print(summary, row.names=FALSE)

write.csv(summary, file.path(OUT, 'ko_summary.csv'), row.names=FALSE)
cat('\nSaved ko_summary.csv\n')
