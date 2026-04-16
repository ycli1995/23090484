
source("/state/partition1/NASDATA1/pipeline/Toolkit_ycli/Eula/v2.0/R/Eula.Seurat.R", chdir = TRUE)

comp = readLines("compare.list")

need.cols = c("log2FC", "P_value", "FDR", "significance")

df = list()
for (i in comp) {
  df[[i]] = ReadTable(paste0(i, ".all.xls")) %>%
    mutate(significance = "nosig") %>%
    mutate(significance = ifelse(log2FC > log(1.5, 2) & FDR < 0.05, "up", significance)) %>%
    mutate(significance = ifelse(log2FC < -log(1.5, 2) & FDR < 0.05, "down", significance)) %>%
    select("ID", starts_with("mean"), !!!need.cols) %>%
    rename(setNames(need.cols, paste(i, need.cols)))
}
df = Reduce(left_join, df) %>%
  relocate(starts_with("mean"), .after = "ID")
colnames(df) = gsub("\\-1", "", colnames(df))

annot = ReadTable("changeid.xls") %>% rename(ID = "Index")

df = left_join(df, annot)
WriteTable(df, "merge_diff.all.xls")

# Venn
library(eulerr)
colors <- c("#F0BB3F", "#B3686F", "#5488AE")
for (i in c("up", "down")) {
  dir.create(i)
  ll = list()
  for (j in comp) {
    ll[[j]] = df[, "ID"][df[, paste(j, "significance")] == i]
  }
  print(str(ll))

  pdf(file = file.path(i, "venn1.pdf"), height=3.5, width=7)
  print(plot(euler(ll), quantities = TRUE, fills = colors, legend = list(side = "right")))
  dev.off()
  pdf(file = file.path(i, "venn2.pdf"), height=3.5, width=7)
  print(plot(venn(ll), quantities = TRUE, fills = colors, legend = list(side = "right")))
  dev.off()
  system(paste("for i in ", file.path(i, "*.pdf"), "; do convert -density 300 $i ${i/pdf/png}; done"))

  venn.sets = get_venn_sets(ll) %>% left_join(df)
  WriteTable(venn.sets, file.path(i, "Venn_sets.xls"))

  inter = venn.sets %>% filter(nchar(set) == max(nchar(set)))
  name = unique(inter$set)
  write(unique(inter$ID), file.path(i, paste0(name, ".glist")))
}


