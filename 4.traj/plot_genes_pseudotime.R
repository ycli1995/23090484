
library(monocle)
library(ComplexHeatmap)
library(dplyr)
library(Seurat)

source("/state/partition1/NASDATA1/pipeline/Toolkit_ycli/Eula/v2.0/R/Eula.Seurat.R", chdir = T)
source("/state/partition1/NASDATA1/pipeline/Toolkit_ycli/Eula/v2.0/R/Eula.monocle2.R", chdir = T)

mnc_obj = readRDS("Trajectory.obj.rds")
obj = Load("obj.Rda")
genes = readLines("genes.list")

branch_point <- mnc_obj@auxOrderingData[[mnc_obj@dim_reduce_type]]$branch_points
branch_name <- FindBranchName(ChangeBranchPointState(mnc_obj), branch_point)

dir.create("1.genes_pseudotime")
for (i in genes) {
  name = FindFeaturesName(obj, i, "name")
  p2 <- plot_genes_branched_pseudotime(mnc_obj[i, ], branch_point = 1, branch_labels = branch_name, color_by = "Samples", min_expr = 0.1)
  colors <- experimentData(mnc_obj)@other$colors[["Samples"]]
  if (length(colors) > 0) {
    p2 <- p2 + scale_color_manual(values = colors)
  }
  ggsave(file.path("1.genes_pseudotime", paste0(name, ".pdf")), p2, width = 5, height = 4.75, limitsize = FALSE)
}

for (b in c(1)) {
  hmp0 = readRDS(paste0("Branch.", b, ".genes_heatmap.Rds"))

  top_annot = hmp0$annotation_col %>% rename("State" = "Cell Type") %>% mutate(State = factor(State))
  left_annot = hmp0$annotation_row %>% arrange(Cluster)
  exp = hmp0$heatmap_matrix[rownames(left_annot), , drop = F]

  gene_names = FindFeaturesName(obj, genes, "name")
  at = match(gene_names, rownames(left_annot))
  labels = gene_names
  label_col = "black"

  #left_cols = list(Cluster = c("#00c4ff", "#ff9289") %>% setNames(levels(left_annot$Cluster)))
  top_cols = list(State = c("grey", "#e36370", "#7591b8") %>% setNames(levels(top_annot$State)))
  cols = hmp0$hmcols

  hmp = Heatmap(
    exp, name = "Expression", col = cols,
    cluster_rows = T,
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    top_annotation = HeatmapAnnotation(df = top_annot, col = top_cols, show_annotation_name = F),
    #left_annot = HeatmapAnnotation(df = left_annot, col = left_cols, show_annotation_name = F, which = "row"),
    right_annot = HeatmapAnnotation(label = anno_mark(at = at, labels = labels, labels_gp = gpar(col = label_col)), show_annotation_name = F, which = "row"),
    row_split = left_annot$Cluster,
    row_title_side = "left",
    row_title_rot = 0,
    row_title = NULL,
    column_split = rep(1:2, each = 100),
    column_title = NULL
  )
  pdf(paste0("Branch.", b, ".heatmap.pdf"), width = 9, height = 10)
  draw(hmp)
  dev.off()

  hmp = Heatmap(
    exp[gene_names, ], name = "Expression", col = cols,
    cluster_rows = T, 
    cluster_columns = F,
    show_row_names = T,
    show_column_names = F,
    top_annotation = HeatmapAnnotation(df = top_annot, col = top_cols, show_annotation_name = F),
    #left_annot = HeatmapAnnotation(df = left_annot, col = left_cols, show_annotation_name = F, which = "row"),
    #right_annot = HeatmapAnnotation(label = anno_mark(at = at, labels = labels, labels_gp = gpar(col = label_col)), show_annotation_name = F, which = "row"),
    #row_split = left_annot$Cluster,
    row_title_side = "left",
    row_title_rot = 0,
    row_title = NULL,
    column_split = rep(1:2, each = 100),
    column_title = NULL
  )
  pdf(paste0("Branch.", b, ".heatmap2.pdf"), width = 9, height = 10)
  draw(hmp)
  dev.off()
}

