
source("/state/partition1/NASDATA1/pipeline/SCPlot/v2.0/bin/enrich_plot_lib.R")
source("/state/partition1/NASDATA1/pipeline/Toolkit_ycli/Eula/v2.0/R/Eula.Seurat.R", chdir = TRUE)

for (i in c("up", "down")) {
  df = ReadTable(file.path(i, "KO.bar_Gradient.xls")) %>% filter(Pvalue < 0.05)
  p = EnrichBarplot(df, top.n = 10)
  h = 10 * 0.2 + 1
  w = 0.08 * max(nchar(df$Descrption)) + 1.5
  ggsave(p, filename = file.path(i, "KO.top10.BarPlot.pdf"), width = w, height = h)

  system(paste("convert -density 300", file.path(i, "KO.top10.BarPlot.pdf"), file.path(i, "KO.top10.BarPlot.png")))
}

