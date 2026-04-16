
source("/state/partition1/NASDATA1/pipeline/Toolkit_ycli/Eula/v2.0/R/Eula.Seurat.R", chdir = T)

all.exps = ReadTable("all_exp.xls", row.names = 1, keep.row.names = F) %>% as.matrix()
metadata = ReadTable("metadata.xls", row.names = 1, keep.row.names = F)[colnames(all.exps), , drop = FALSE]
changeid = ReadTable("changeid.xls")

use.spots = Matrix::colSums(all.exps) > 0
all.exps = all.exps[, use.spots]
metadata = metadata[use.spots, , drop = FALSE]

all.diff = ReadTable("diff.list", header = FALSE) %>% as.matrix()

for (i in seq_len(nrow(all.diff))) {
  all.groups = all.diff[i, , drop = TRUE]
  compare = paste(all.groups, collapse = "-vs-")

  cells1 = rownames(metadata)[metadata$Sample == all.groups[1]]
  cells2 = rownames(metadata)[metadata$Sample == all.groups[2]]

  exps = all.exps[, c(cells1, cells2), drop = FALSE]
  groups = c(rep("group1", length(cells1)), rep("group2", length(cells2)))

  pseudo.count = 0.000001
  mean1 = Matrix::rowMeans(all.exps[, cells1])
  mean2 = Matrix::rowMeans(all.exps[, cells2])
  mean1 = ifelse(mean1 == 0, pseudo.count, mean1)
  mean2 = ifelse(mean2 == 0, pseudo.count, mean2)
  log2fc = log2(mean2 / mean1)

  pval = vapply(seq_len(nrow(all.exps)), function(ii) {
    if (ii %% 100 == 0) print(ii)
    wilcox.test(exps[ii, ] ~ groups)$p.value
  }, numeric(1))
  fdr = p.adjust(pval, method = "fdr")

  out = data.frame(
    ID = rownames(exps),
    mean1 = mean1, 
    mean2 = mean2,
    log2FC = log2fc,
    P_value = pval,
    FDR = fdr
  )
  colnames(out)[2:3] = paste0("mean_", all.groups)

  out = out %>% left_join(changeid %>% rename(ID = "Index"))
  WriteTable(out, paste0(compare, ".all.annot.xls"))

  filtered = out %>% filter(abs(log2FC) > log2(1.5) & FDR < 0.05)
  WriteTable(filtered, paste0(compare, ".filtered.annot.xls"))

  dir.create("enrich")
  glist = filtered[, c("ID", "log2FC")]
  WriteTable(glist, paste0("enrich/", compare, ".glist"), col.names = F)
}

