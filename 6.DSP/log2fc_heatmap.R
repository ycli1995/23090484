
source("/state/partition1/NASDATA1/pipeline/Toolkit_ycli/Eula/v2.0/R/Eula.Seurat.R", chdir = TRUE)

df1 = ReadTable("DSP.xls") %>%
  filter(nchar(set) == max(nchar(set)))
df2 = ReadTable("scRNA.xls") %>%
  filter(nchar(set) == max(nchar(set))) %>%
  select(!Cluster)

genes = readLines("venn/DSP__scRNA.glist")

m1 = df1 %>% 
  filter(GeneName %in% genes) %>%
  tibble::column_to_rownames("GeneName") %>%
  select(matches("log2FC")) %>%
  as.matrix()
colnames(m1) = gsub(".*-vs-", "", gsub(" .*", "", colnames(m1)))

m2 = df2 %>% 
  filter(GeneName %in% genes) %>%
  tibble::column_to_rownames("GeneName") %>%
  select(matches("log2FC")) %>%
  as.matrix()
colnames(m2) = gsub(".*-vs-", "", gsub(" .*", "", colnames(m2)))

m1 = abs(m1[genes, , drop = FALSE])
m2 = abs(m2[genes, , drop = FALSE])
m1 = apply(m1, 1, minMax) %>% t()
m2 = apply(m2, 1, minMax) %>% t()

library(circlize)
library(ComplexHeatmap)

ht2 = Heatmap(m2, cluster_columns = F, cluster_rows = T)
pdf(tempfile())
ht2 = draw(ht2)
dev.off()
order = row_order(ht2)

genes = rownames(m2)[order]

#m1 = abs(m1[genes, , drop = FALSE])
#m2 = abs(m2[genes, , drop = FALSE])

#m1 = minMax(m1)
#m2 = minMax(m2)
#m1 = t(scale(t(m1), center = F))
#m2 = t(scale(t(m2), center = F))
#m1[m1 > 5] = 5
#m2[m2 > 5] = 5

mats = list(DSP = m1, scRNA = m2)
colors = list(
  DSP = circlize::colorRamp2(c(0, max(m1)), c("white", "#d05f27")),
  scRNA = circlize::colorRamp2(c(0, max(m2)), c("white", "#407fc2"))
)

size = 0.12
ht = NULL
for (i in names(mats)) {
  m = mats[[i]]
  if (nrow(m) > 30) size = 0.07
  ht = ht + Heatmap(
    m, 
    col = colors[[i]], 
    name = paste(i, "log2(fc)"), 
    column_title = i, 
    show_row_names = nrow(m) < 30,
    cluster_rows = F, 
    cluster_columns = F,
    heatmap_legend_param = list(labels_gp = gpar(fontsize = 0), direction = "horizontal", legend_width = unit(20, "mm"))
  )
} 
h = max(4, nrow(m) * size + 1.75)
w = 6
outfile = paste0("log2FC", ".heatmap.pdf")
pdf(outfile, height = h, width = w)
draw(ht, heatmap_legend_side = "top")
dev.off()
system(paste("convert -density 300", outfile, gsub("pdf", "png", outfile)))

