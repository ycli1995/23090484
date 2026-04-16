
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

source("/public2/Bio/pipeline/SingleCell_Collections/SCellWare/v3.1/R/Seurat_lib.R", chdir = T)

old_env <- environment(.PlotFeaturePlot)
old_env$.PlotFeaturePlot = function(object, features, outfile = NULL, reduction = NULL, is.use.name = TRUE, color.high = "blue", color.low = "lightgrey", show.cluster.label = FALSE, nCol = NULL, plot.basic.size = 4, group.by = "seurat_clusters", is.combine = TRUE, cols = NULL, ... ) {
		if ( show.cluster.label ) Idents(object) <- group.by
		if ( is.null(cols) ) cols <- c(color.low, color.high)
		plots <- FeaturePlot(object, features = features, order = TRUE, reduction = reduction, combine = FALSE, label = show.cluster.label, cols = cols, ...)
		for (i in seq_along(plots)) {
			plots[[i]] = plots[[i]] + scale_color_gradientn(limits = c(0, 4), colors = cols)
		}
		if ( is.use.name ) {
				name <- FindFeaturesName(object, features)
				for ( i in seq(plots) ){
						plots[[i]] <- plots[[i]] + ggtitle(name[i])
				}
				names(plots) <- name
		} else {
				names(plots) <- features
		}
		if ( is.null(nCol) ) nCol <- ceiling(sqrt(length(features)))
		nRow <- ceiling(length(features) / nCol)

		p <- wrap_plots(plots, ncol = nCol) & dot_theme_default()

		if ( is.null(outfile) ) {
				return(p)
		} else {
				ggsave(p, file = outfile, width = plot.basic.size * (6/5) * nCol, height = plot.basic.size * nRow, limitsize = FALSE)
		}
}
detach(name = lib_conf$package.name, character.only = TRUE)
attach(old_env, name = lib_conf$package.name)

obj = Load("obj_renamed.Rda")

genes = readLines("genes.list")

## umap
dir.create("umap")
#dir.create("tsne")

cols = c("#B1DAE9", "#FE0F06")
PlotFeaturePlot(obj, FindFeaturesID(obj, genes), reduction = "umap", is.combine = FALSE, outpref = "umap/All", cols = cols)
#PlotFeaturePlot(obj, FindFeaturesID(obj, genes), reduction = "tsne", is.combine = FALSE, outpref = "tsne/All", cols = cols)
#for (i in levels(obj$Groups)) {
#	obj2 = subset(obj, Groups %in% i)
#	PlotFeaturePlot(obj2, FindFeaturesID(obj, genes), reduction = "umap", is.combine = FALSE, outpref = paste0("umap/", i), cols = cols)
#	PlotFeaturePlot(obj2, FindFeaturesID(obj, genes), reduction = "tsne", is.combine = FALSE, outpref = paste0("tsne/", i), cols = cols)
#}

p = PlotDotPlot(obj, FindFeaturesID(obj, genes)) + 
  scale_colour_gradientn(colors = colorRampPalette(c("#B1DAE9", "#FE0F06"))(100)) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 0)
  ) + 
  theme(legend.box = "horizontal")
p$data$id = factor(as.character(p$data$id), levels = rev(levels(p$data$id))) 
p$data$features.plot = factor(as.character(p$data$features.plot), levels = rev(levels(p$data$features.plot)))
levels(p$data$features.plot) = gsub("\\s\\(.*", "", levels(p$data$features.plot))

w <- max(7, ceiling(length(genes)) * 0.35 + 2)
h = max(2.75, length(unique(obj@meta.data[["seurat_clusters"]])) * 0.4)
ggsave(p, file = "DotPlot.pdf", width = w, height = h, limitsize = FALSE )

p = VlnPlot(obj, features = FindFeaturesID(obj, genes), group.by = "seurat_clusters", pt.size = 0, stack = TRUE, cols = obj@misc$color.cluster, fill.by = "ident", flip = FALSE)
levels(p$data$feature) = FindFeaturesName(obj, levels(p$data$feature), col = "name")
p$data$ident = factor(as.character(p$data$ident), rev(levels(p$data$ident)))
p = p +
  labs(x = "", y = "") +
  NoLegend() +
  theme(axis.ticks.x.bottom = element_blank(), axis.text.x.bottom = element_blank(), strip.text.x = element_text(hjust = 1, face = "plain"))

w = 5.5
h = 0.4 * length(genes) + 1.25
ggsave(p, filename = "Stack_VlnPlot.pdf", width = h, height = w)

CalAvgExp(obj)
obj$cluster_group = paste0(obj$seurat_clusters, "_", obj$Groups)
obj$cluster_group = factor(obj$cluster_group, levels = as.vector(sapply(levels(obj$seurat_clusters), paste0, "_", levels(obj$Groups))))
CalAvgExp(obj, outfile = "AllGene.Group.avg_exp.xls", group.by = "cluster_group")

write(FindFeaturesID(obj, genes), "genes.list")

