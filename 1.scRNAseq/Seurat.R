
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(showtext)
showtext_auto(enable = TRUE)
library(RColorBrewer)
source("Seurat_lib.R", chdir = T)

### Let's shake it
obj = Load("obj_renamed.Rda")

source("clusterCornerAxes.R")

plotCluster = function(obj, reduction = "umap", p1 = "orig.ident", p2 = "seurat_clusters", outpfx = "UMAP", ...) {
  plist = list()

  g1 = obj[[]][[p1]]
  plist[[1]] = clusterCornerAxes(obj, reduction = reduction, clusterCol = p1, pSize = 0.01) +
    labs(color = "") + 
    scale_color_manual(values = obj@misc$colors[[p1]][levels(g1)])
  char_len1 = max(nchar(as.character(g1))) * 0.05 + 0.05

  g2 = obj[[]][[p2]]
  plist[[2]] = clusterCornerAxes(obj, reduction = reduction, clusterCol = p2, pSize = 0.01) +
    labs(color = "") + 
    scale_color_manual(values = obj@misc$colors[[p2]][levels(g2)])
  char_len2 = max(nchar(as.character(g2))) * 0.05 + 0.05
  p = wrap_plots(plist, nrow = 1)
  h = 6
  w = 6 + char_len1 + 6 + char_len2
  ggsave(p, filename = paste0(outpfx, ".pdf"), height = h, width = w, limitsize = FALSE)
}

plotCluster(obj)
plotCluster(obj, p1 = "Groups", outpfx = "UMAP.Groups")

for (i in levels(obj$orig.ident)) {
  obj2 = subset(obj, orig.ident %in% i)
  obj2@meta.data = droplevels(obj2@meta.data)
  i2 = gsub("\\/|\\s", "_", i)
  plotCluster(obj2, outpfx = paste0("UMAP.", i2))
}

for (i in levels(obj$Groups)) {
  obj2 = subset(obj, Groups %in% i)
  obj2@meta.data = droplevels(obj2@meta.data)
  i2 = gsub("\\/|\\s", "_", i)
  plotCluster(obj2, p1 = "Groups", outpfx = paste0("UMAP.Groups.", i2))
}

## Pie plot
cell = data.frame(
	Cells = Cells(obj),
	Cluster = obj$seurat_clusters,
	Groups = obj$Groups
)
color = obj@misc$color.cluster
library(ggrepel)
for (i in levels(cell$Groups)) {
  df = cell[cell$Groups == i, ]
  pdata = df %>% 
    group_by(Cluster) %>% 
    summarise(Freq = length(Cells), .groups = 'drop') %>%
    mutate(fraction = Freq / sum(Freq)) %>%
    mutate(ymax = cumsum(fraction)) %>%
    mutate(ymin = c(0, ymax[-nrow(.)])) %>%
    mutate(Cluster2 = paste0(round(fraction * 100, digits = 2), "%"))

  p = ggplot(pdata, aes(fill = Cluster, ymax = ymax, ymin = ymin, xmax = 3.5, xmin = 2, col = Cluster)) +
    geom_rect(color = "white") + coord_polar('y') + xlim(c(0, 4)) +
    scale_fill_manual(values = color[pdata$Cluster]) + 
#    scale_color_manual(values = color[pdata$Cluster]) +
#    geom_label_repel(aes(label = Cluster2, x = 3.5, y = (ymax + ymin)/2), show_guide  = FALSE, color = 'black', size = 4,  max.overlaps = Inf) +
    labs(fill = 'Cluster', color = 'Cluster') +
    theme_void() +
    theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  i2 = gsub("\\s|\\/", "_", i)
  ggsave(paste0(i2, '.cluster_pie.pdf'), p, w = 7, h = 5)
}

## Bar plot
.PlotClusterStat = function(object, stat.what = "seurat_clusters", group.by = "orig.ident", color.st = NULL, color.gb = NULL, outpref = NULL, ...){
		metadata <- object@meta.data
		if ( is.null(color.st) ) {
				if ( "misc" %in% slotNames(object) && exists(stat.what, object@misc) ) {
						color.st <- object@misc[[stat.what]]
				} else {
						color.st <- switch(stat.what,
								"seurat_clusters" = object@misc$color.cluster,
								"orig.ident" = object@misc$color.sample,
								"Groups" = object@misc$color.group
						)
				}
		}
		if ( is.null(color.gb) ) {
				if ( "misc" %in% slotNames(object) && exists(group.by, object@misc) ) {
						color.gb <- object@misc[[group.by]]
				} else {
						color.gb <- switch(group.by, 
								"seurat_clusters" = object@misc$color.cluster,
								"orig.ident" = object@misc$color.sample,
								"Groups" = object@misc$color.group
						)
				}
		}
		name.st <- switch(stat.what, "seurat_clusters" = "Cluster", "orig.ident" = "Samples", stat.what)
		name.gb <- switch(group.by,  "seurat_clusters" = "Cluster", "orig.ident" = "Samples", group.by)
		stat.what <- as.name(stat.what)
		group.by  <- as.name(group.by)
		stat_sample <- metadata %>%
				group_by(!! name.gb := !! group.by, !! name.st := !! stat.what) %>%
				summarise("Number of cells" = n())
		p <- list()
		p[["by"]] <- ggplot(stat_sample, aes_(x = as.name(name.gb), y = ~ `Number of cells`, fill = as.name(name.st)))
		p[["in"]] <- ggplot(stat_sample, aes_(x = as.name(name.st), y = ~ `Number of cells`, fill = as.name(name.gb)))
		if ( ! is.null(color.st) ) p[["by"]] <- p[["by"]] + scale_fill_manual(values = color.st)
		if ( ! is.null(color.gb) ) p[["in"]] <- p[["in"]] + scale_fill_manual(values = color.gb)
		geom_stack <- geom_bar(stat = "identity", position = 'stack')
		geom_fill  <- geom_bar(stat = "identity", position = "fill" )
		if ( is.null(outpref) ) {
				outpref <- paste0( name.st, ".stat")
		}
		for ( i in names(p) ) {
				p[[i]] <- p[[i]] + guides(fill=guide_legend(ncol = 2)) + bar_theme_default() + theme(plot.margin = margin(l = 20))
				if (i == "by") h = length(unique(stat_sample[, name.gb])) * 1 + 4
				if (i == "in") h = length(unique(stat_sample[, name.st])) * 1 + 4
				w = 5 + max(nchar(unique(stat_sample[, name.gb]))) * 0.055 + max(nchar(unique(stat_sample[, name.st]))) * 0.055
				ggsave( p[[i]] + geom_stack + coord_flip(), file = paste0( outpref, ".", i, name.gb, ".pdf"), height = h, width = w)
				ggsave( p[[i]] + geom_fill + ylab("Fraction of Cells") + coord_flip(),  file = paste0( outpref, ".", i, name.gb, ".pct.pdf"), height = h, width = w )
		}
}
StatCluster = function(object, group.by = "orig.ident", outpref = "Cluster.stat", stat.what = "seurat_clusters", assay = DefaultAssay(object), suffix = 'Cell', ...){
	.StatCluster(object, stat.what = stat.what, outpref = outpref, assay = assay, suffix = suffix)
	.StatCluster_by(object, stat.what = stat.what, group.by = group.by, outpref = outpref)
	.PlotClusterStat(object, stat.what = stat.what, group.by = group.by, outpref = outpref, ...)
}
StatCluster(obj)
StatCluster(obj, group.by = "Groups", outpref = "Cluster.stat")

