cc <- readRDS('cellchat_output_subtype_seed123/cellchat_Strong.rds')
cat('Idents:', levels(cc@idents), '\n')
cat('Treg in idents:', any(c('Treg_FOXP3','Treg_CCR8','Treg_MKI67') %in% levels(cc@idents)), '\n')
cat('Net prob dims:', dim(cc@net$prob), '\n')
cat('Number of LRs:', length(cc@net$LRs), '\n')
cat('First few LRs:', head(cc@net$LRs, 5), '\n')

# Check if any Treg source has non-zero prob
treg_groups <- c('Treg_FOXP3','Treg_CCR8','Treg_MKI67')
idents <- levels(cc@idents)
treg_idx <- which(idents %in% treg_groups)
cat('Treg indices:', treg_idx, '\n')
cat('Total ident levels:', length(idents), '\n')

found <- 0
for (si in treg_idx) {
  for (ti in 1:length(idents)) {
    s <- sum(cc@net$prob[si, ti, , drop = FALSE])
    if (s > 0) {
      cat('  Non-zero:', idents[si], '->', idents[ti], ':', s, '\n')
      found <- found + 1
      if (found >= 10) break
    }
  }
  if (found >= 10) break
}
if (found == 0) cat('No non-zero Treg source interactions found in Strong\n')

# Also check Weak
cc2 <- readRDS('cellchat_output_subtype_seed123/cellchat_Weak.rds')
idents2 <- levels(cc2@idents)
treg_idx2 <- which(idents2 %in% treg_groups)
found2 <- 0
for (si in treg_idx2) {
  for (ti in 1:length(idents2)) {
    s <- sum(cc2@net$prob[si, ti, , drop = FALSE])
    if (s > 0) {
      cat('  Weak Non-zero:', idents2[si], '->', idents2[ti], ':', s, '\n')
      found2 <- found2 + 1
      if (found2 >= 10) break
    }
  }
  if (found2 >= 10) break
}
if (found2 == 0) cat('No non-zero Treg source interactions found in Weak\n')
