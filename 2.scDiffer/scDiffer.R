#args <- commandArgs()
#bin  <- dirname(normalizePath(sub('--file=', '',  args[grep('--file=', args)])))
#args <- args[-seq(grep("--args", args))]
args <- commandArgs(T)

file   <- args[1]
outdir  <- args[2]
add_lib <- args[3]
lib.loc <- args[4]

if ( ! is.na(lib.loc) ) {
		.libPaths(lib.loc)
		options(stringsAsFactors = TRUE)
}

parameter <- yaml::yaml.load_file(file, handlers = list("bool#no" = function(x){if ( x %in% c("false", "FALSE") ) FALSE else x}, "bool#yes" = function(x){if ( x %in% c("true", "TRUE") ) TRUE else x}))

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

#library(future)
#options(future.globals.maxSize = 100 * 1024^3)
#plan("multicore", workers = 4)
#plan("sequential")

source(add_lib, chdir = T)


setwd(outdir)

obj <- Load(parameter$data$object)

if ( ! is.null(parameter$data$assay) ) {
		DefaultAssay(obj) <- parameter$data$assay
}

if ( !is.null(parameter$cluster_col) ) {
		cluster_col <- parameter$cluster_col
} else {
		cluster_col <- 'seurat_clusters'
}

if ( !is.null(parameter$rename_clusters) ) {
		new_cluster <- unlist.rev(parameter$rename_clusters)
		old_cluster <- `names<-`(levels(obj@meta.data[[cluster_col]]), levels(obj@meta.data[[cluster_col]]))
		new_cluster <- c(new_cluster, old_cluster[!names(old_cluster) %in% names(new_cluster)])
		levels(obj@meta.data[[cluster_col]]) <- new_cluster[levels(obj@meta.data[[cluster_col]])]
		Idents(obj) <- cluster_col
		obj@misc$color.cluster <- SetColor(Idents(obj))
} else {
		Idents(obj) <- cluster_col
}

if ( !is.null(parameter$cluster.use) ) {
		obj <- SubsetObj(obj, cluster = parameter$cluster.use, cluster.name = cluster_col)
}
if ( !is.null(parameter$cell.use) ) {
		cells <- readLines(parameter$cell.use)
		obj <- SubsetObj(obj, cells = cells)
}
obj@meta.data <- droplevels(obj@meta.data)

if ( !is.null(parameter$feature.use) ) {
		features <- readLines(parameter$feature.use)
		obj <- obj[features,]
}

group.data <- FindGroupIndex(obj, parameter)
CalAvgByGrpInCls(obj, group.data, group_names = unique(unlist(parameter$differ)), diff_type = parameter$diff_type)

marker <- DoFindGroupMarkers(obj, group.data, parameter)

save(marker, file = "marker.Rda")
WriteDifferMarker(marker, "all", object = obj)
StatDeGene(marker, group.by = "contrast")
PlotDeGeneVolcano(marker, parameter)

diff.marker <- lapply(marker, function(x) filter(x, sig != "nosig"))
WriteDifferMarker(diff.marker, "diff", object = obj)

PlotContrastHeatmap(obj, diff.marker, group.data, top.num = parameter$plots$top, cluster.name = cluster_col)
PlotContrastDotPlot(obj, diff.marker, group.data, top.num = parameter$plots$top, cluster.name = cluster_col)
PlotContrastFeature(obj, diff.marker, group.data, top.num = parameter$plots$top, reduction = 'tsne', cluster.name = cluster_col)
PlotContrastFeature(obj, diff.marker, group.data, top.num = parameter$plots$top, reduction = 'tsne', is.demo = T, cluster.name = cluster_col)
PlotContrastFeature(obj, diff.marker, group.data, top.num = parameter$plots$top, reduction = 'umap', cluster.name = cluster_col)
PlotContrastFeature(obj, diff.marker, group.data, top.num = parameter$plots$top, reduction = 'umap', is.demo = T, cluster.name = cluster_col)

