FindGroupIndex <- function(object, parameter = list(), slot = "counts", assay = NULL,
				sample_col = IfNull(parameter$sample_col, "orig.ident"), cluster_col = IfNull(parameter$cluster_col, "seurat_clusters"),
				by_sample = parameter$group$by_sample,
				by_cluster = parameter$group$by_cluster,
				by_gene = parameter$group$by_gene
				){
		metadata <- object@meta.data %>% select(sample = !! sample_col)
		for ( i in unique(metadata$sample)) {
				metadata <- metadata %>% mutate(!! i := sample %in% i )
		}
		for ( i in names(by_sample) ){
				metadata <- metadata %>% mutate(!! i := sample %in% by_sample[[i]])
		}
		metadata <- metadata %>% select(- sample)

		if ( ! is.null(by_cluster) ) {
				metadata2 <- object@meta.data %>% select(cluster = !! cluster_col)
				for ( i in names(by_cluster) ){
						metadata2 <- metadata2 %>% mutate(!! i := cluster %in% by_cluster[[i]])
				}
				metadata2 <- metadata2 %>% select(- cluster)
				metadata <- cbind(metadata, metadata2)
		}

		for ( i in names(by_gene) ) {
				gene_id <- by_gene[[i]][[1]]
				gene_id <- FindFeaturesID(object, gene_id)
				gene_thres <- by_gene[[i]][[2]]
				data <- GetAssayData(object, slot = slot, assay = assay)
				metadata <- metadata %>% mutate(!! i := data[gene_id, ] >= gene_thres[[1]] & data[gene_id, ] <= gene_thres[[2]])
		}

		rownames(metadata) <- colnames(object)
		return(metadata)
}

FindGroupMarkers <- function (object, group.data, contrast = NULL,
		return.thresh = 0.01, use.qvalue = FALSE,
		assay = NULL, features = NULL, logfc.threshold = 0.25, base = 2,
		test.use = "MAST", slot = "data", min.pct = 0.1, min.diff.pct = -Inf, 
		verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
		random.seed = 1, latent.vars = NULL, min.cells.feature = 3, 
		min.cells.group = 3, pseudocount.use = 1, is.bulk = TRUE, is.cluster = TRUE, ...) {
		
		if ((test.use == "roc") && (return.thresh == 0.01)) {
				return.thresh <- 0.7
		}

		cells.1 <- rownames(group.data)[group.data[[contrast[[1]]]]]
		cells.2 <- rownames(group.data)[group.data[[contrast[[2]]]]]
		if ( length(cells.1) == 0 || length(cells.2) == 0 ) {
				message( "===> cell num : \n",
						"        ", contrast[[1]], " : ", length(cells.1), "\n",
						"        ", contrast[[2]], " : ", length(cells.2), "\n",
						"        [SKIP]")
				return(NULL)
		}

		data.1 <- CalAvgExp(object[, cells.1], is.return = TRUE, is.bulk = is.bulk) %>%
				as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>%
				reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "mean.1") %>%
				mutate(cluster = as.factor(cluster))
		data.2 <- CalAvgExp(object[, cells.2], is.return = TRUE, is.bulk = is.bulk) %>%
				as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>% 
				reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "mean.2") %>%
				mutate(cluster = as.factor(cluster))
		pct.1 <- CalPctExp(object[, cells.1], is.return = TRUE, is.bulk = is.bulk) %>% round(digits = 3) %>%
				as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>% 
				reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "pct.1") %>%
				mutate(cluster = as.factor(cluster))
		pct.2 <- CalPctExp(object[, cells.2], is.return = TRUE, is.bulk = is.bulk) %>% round(digits = 3) %>%
				as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>% 
				reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "pct.2") %>%
				mutate(cluster = as.factor(cluster))
		c1 <- table(Idents(object[,cells.1]))
		c2 <- table(Idents(object[,cells.2]))
		if ( is.bulk ) {
				c1 <- c(bulk = length(cells.1), c1)
				c2 <- c(bulk = length(cells.2), c2)
		}
		data <- full_join(x = data.1, y = data.2, by = c("gene", "cluster")) %>%
				left_join(y = pct.1, by = c("gene", "cluster") ) %>%
				left_join(y = pct.2, by = c("gene", "cluster") ) %>%
				mutate(cells.1 = c1[cluster], cells.2 = c2[cluster]) %>%
				mutate(cluster = factor(cluster, levels = if ( is.bulk ) c("bulk", levels(object)) else levels(object))) %>%
				mutate(avg_logFC = log2((mean.2 + pseudocount.use)/(mean.1 + pseudocount.use)) / log2(base)) %>%
#				mutate(avg_logFC = log2(mean.2/mean.1) / log2(base)) %>%
				mutate(p_val = 1, p_val_adj = 1, sig = "nosig")
				
		data[is.na(data)] <- 0
		rm(data.1, data.2, pct.1, pct.2, c1, c2)

#		idents.all <- if ( is.bulk ) c("bulk", levels(droplevels(object))) else levels(droplevels(object))
		if ( is.bulk && is.cluster ) {
				idents.all <- c("bulk", levels(droplevels(object)))
		} else if ( is.bulk ) {
				idents.all <- "bulk"
		} else if ( is.cluster ) {
				idents.all <- levels(droplevels(object))
		}
		data <- data %>% filter(cluster %in% idents.all)

		if ( is.null(features) ) {
				features.use <- data %>% filter( #(mean.1 + pseudocount.use) * (mean.2 + pseudocount.use) != 0 &
												mean.1 + mean.2 > 0 &
												(pct.1 > min.pct | pct.2 > min.pct) &
												abs(avg_logFC) > logfc.threshold / log(base))
				features.use <- split(features.use, features.use$cluster)
				features.use <- lapply(features.use, function(x) as.character(x$gene))
				logfc.threshold <- 0
#				pseudocount.use <- 0
		}else{
				features.use <- lapply(idents.all, function(x) features)
		}

		genes.de <- list()
		messages <- list()
		for (i in seq(idents.all)) {
				if (verbose) message("Calculating cluster ", idents.all[i])


				cells.idents <- if ( idents.all[[i]] == "bulk" ){ 
						WhichCells(object)
				} else {
						WhichCells(object, idents = idents.all[[i]])
				}
				cells.1.tmp <- intersect(cells.2, cells.idents) # switch order from T-vs-C to C-vs-T
				cells.2.tmp <- intersect(cells.1, cells.idents) # switch order from T-vs-C to C-vs-T
				if ( packageVersion("Seurat") >= as.numeric_version("4.0.0") ) {
						genes.de[[i]] <- tryCatch(expr = {
								FindMarkers(object = object, ident.1 = cells.1.tmp, ident.2 = cells.2.tmp,
											assay = assay,
											features = features.use[[idents.all[[i]]]], logfc.threshold = logfc.threshold, 
											test.use = test.use, slot = slot, min.pct = min.pct, 
											min.diff.pct = min.diff.pct, verbose = verbose, 
											only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
											random.seed = random.seed, latent.vars = latent.vars, 
											min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
											pseudocount.use = pseudocount.use, 
											base = exp(1), ## keep same as Seurat < 4.0 
											...)
								}, error = function(cond) { return(cond$message) })
				} else {
						genes.de[[i]] <- tryCatch(expr = {
								FindMarkers(object = object, ident.1 = cells.1.tmp, ident.2 = cells.2.tmp,
											assay = assay,
											features = features.use[[idents.all[[i]]]], logfc.threshold = logfc.threshold, 
											test.use = test.use, slot = slot, min.pct = min.pct, 
											min.diff.pct = min.diff.pct, verbose = verbose, 
											only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
											random.seed = random.seed, latent.vars = latent.vars, 
											min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
											pseudocount.use = pseudocount.use, 
											...)
								}, error = function(cond) { return(cond$message) })
				}
				if (class(x = genes.de[[i]]) == "character") {
						messages[[i]] <- genes.de[[i]]
						genes.de[[i]] <- NULL
				}
		}

		gde.all <- data.frame()
		for (i in seq(idents.all)) {
				if (is.null(x = unlist(x = genes.de[i]))) {
						next
				}
				gde <- genes.de[[i]]
				if (nrow(x = gde) > 0) {
						if (test.use == "roc") {
								gde <- subset(x = gde, subset = (myAUC > return.thresh | myAUC < (1 - return.thresh)))
#						} else if (is.null(x = node) || test.use %in% c("bimod", "t")) {
						} else {
								gde <- gde[order(gde$p_val_adj, gde$p_val, -gde[, 2]), ]
								if ( use.qvalue ) {
										gde <- subset(x = gde, subset = p_val_adj < return.thresh)
								} else {
										gde <- subset(x = gde, subset = p_val < return.thresh)
								}
						}
				
#						if (nrow(x = gde) > 0) {
								genes.de[[i]]$sig <- "nosig"
								genes.de[[i]][rownames(gde)[gde$avg_logFC > 0], "sig"] <- "up"
								genes.de[[i]][rownames(gde)[gde$avg_logFC < 0], "sig"] <- "down"	
								genes.de[[i]]$cluster <- idents.all[i]
								genes.de[[i]]$gene <- rownames(x = genes.de[[i]])
								gde.all <- rbind(gde.all, genes.de[[i]])
#						}
				}
		}
		if ( nrow(x = gde.all) == 0) {
				warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
				data <- data %>%
						mutate(sig = factor(sig, levels = c("up", "nosig", "down")), 
						       cluster = factor(cluster, levels = levels(data$cluster)), 
							   contrast = factor(paste0(contrast, collapse = "-vs-")) ) %>%	
						select(cluster, gene, everything())
				return(data)
		}
		gde.all <- gde.all %>%
				select(cluster, gene, p_val, p_val_adj, sig) %>%
				full_join(y = data, by = c("gene", "cluster")) %>%
				mutate( p_val     = if_else(is.na(p_val.x),     p_val.y,     p_val.x),
						p_val_adj = if_else(is.na(p_val_adj.x), p_val_adj.y, p_val_adj.x),
						sig       = if_else(is.na(sig.x),       sig.y,       sig.x)) %>%						
				mutate( sig     = factor(sig, levels = c("up", "nosig", "down")),
						cluster = factor(cluster, levels = levels(data$cluster)),
						contrast  = factor(paste0(contrast, collapse = "-vs-")) ) %>%
				select(-p_val.x, -p_val.y, -p_val_adj.x, -p_val_adj.y, -sig.x, -sig.y)

		if ((only.pos) && nrow(x = gde.all) > 0) {
				return(subset(x = gde.all, subset = gde.all[, "avg_logFC"] > 0))
		}
		if (nrow(x = gde.all) == 0) {
				warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
		}
		if (length(x = messages) > 0) {
				warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
				for (i in 1:length(x = messages)) {
						if (!is.null(x = messages[[i]])) {
								warning("When testing cluster ", idents.all[i], " :\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
						}
				}
		}
		return(gde.all)
}

DoFindGroupMarkers <- function(object, group.data, parameter = list(),
				differ_groups   = parameter$differ,
				return.thresh   = IfNull(eval(parse(text = parameter$FindMarkers$pvalue)), 0.01),
				use.qvalue      = IfNull(parameter$FindMarkers$use.q, FALSE),
				logfc.threshold = IfNull(eval(parse(text = parameter$FindMarkers$log2fc)) / log2(exp(1)), 0.25),
				test.use        = IfNull(parameter$FindMarkers$method, "MAST"),
				min.pct         = IfNull(eval(parse(text = parameter$FindMarkers$min_pct)), 0.1),
				pseudocount.use = IfNull(eval(parse(text = parameter$FindMarkers$pseudocount.use)), 1),
				diff_type       = IfNull(parameter$diff_type, "all"),
				...){
		if ( diff_type == "all" )  {
				is.bulk <- TRUE
				is.cluster <- TRUE
		} else if ( diff_type == "no_bulk" ) {
				is.bulk <- FALSE
				is.cluster <- TRUE
		} else if ( diff_type == "only_bulk" ) {
				is.bulk <- TRUE
				is.cluster <- FALSE
		} else {
				stop()
		}
		if ( future::nbrOfWorkers() > 1 ) {
				marker <- future.apply::future_lapply(differ_groups, function(i)
								FindGroupMarkers(object, group.data, contrast = i,
										return.thresh = return.thresh,
										use.qvalue = use.qvalue,
										logfc.threshold = logfc.threshold,
										test.use = test.use,
										min.pct = min.pct, 
										pseudocount.use = pseudocount.use,
										is.bulk = is.bulk, is.cluster = is.cluster,
										...
										)
						)
				names(marker) <- sapply(differ_groups, paste, collapse = "-vs-")
				index <- sapply(marker, is.null)
				marker <- marker[!index]
		}else{
				marker <- list()
				for ( i in differ_groups ) {
						name <- paste(i, collapse = "-vs-")
						marker[[name]] <- FindGroupMarkers(object, group.data, contrast = i,
										return.thresh = return.thresh,
										use.qvalue = use.qvalue,
										logfc.threshold = logfc.threshold,
										test.use = test.use,
										min.pct = min.pct, 
										pseudocount.use = pseudocount.use,
										is.bulk = is.bulk, is.cluster = is.cluster,
										...
										)
				}
		}
		return(marker)
}


WriteDifferMarker <- function(marker, add_name = NULL, object = NULL) {
		for ( contrast_name in names(marker) ) {
				contrast <- strsplit(contrast_name, "-vs-")[[1]]
				if ( ! is.null(object) ) {
						marker[[contrast_name]]$gene <- ChangeOUTName(marker[[contrast_name]]$gene, object@misc$fdata)
						marker[[contrast_name]]$name <- FindFeaturesName(object, marker[[contrast_name]]$gene, "name")
				} else {
						marker[[contrast_name]]$name <- marker[[contrast_name]]$gene
				}
				data <- marker[[contrast_name]] %>%
						select(Cluster = cluster, GeneID = gene, GeneName = name,
								!! contrast[[1]] := mean.1, !! contrast[[2]] := mean.2,
								cells.1, cells.2, pct.1, pct.2,
								log2FC = avg_logFC, p_value = p_val, p_val_adjust = p_val_adj, significance = sig)
				WriteTable(data, file = paste(c("DifferMarker", contrast_name, add_name, "xls"), collapse = "." ))
		}
		all.list <- do.call(rbind, lapply(names(marker), function(i) { if ( nrow(marker[[i]]) == 0 ) return(NULL); cbind(contrast = i, marker[[i]][,c("cluster", "gene")]) } ) )
		WriteTable(all.list, file = "all.differ.list")
}

StatDeGene <- function(marker, group.by = c("contrast", "cluster"), do.plot = TRUE){
		group.by <- match.arg(group.by)
		if ( is.null(marker) || length(marker) == 0 ) return(NULL)
		stat <- do.call(rbind, lapply(marker, function(x) x %>%
								group_by(cluster, contrast, cells.1, cells.2, sig) %>%
								summarise(count = n()) %>%
								as.data.frame() ) )
#		stat2 <- reshape2::dcast(stat, cluster + contrast ~ sig, value.var = "count", fill = 0, drop = F) %>% mutate(total = up + down)
		stat2 <- reshape2::dcast(stat, ... ~ sig, value.var = "count", fill = 0, drop = T) %>%
				mutate(up = if(exists("up")) up else 0, down = if(exists("down")) down else 0) %>% 
				select(-nosig) %>% mutate(total = up + down) %>% filter(cells.1 * cells.2 != 0)
		for ( i in levels(stat2[[group.by]]) ) {
				name <- i
				if ( group.by == "contrast" ) {
#						name <- paste0("Contrast.", name)
				} else if ( group.by == "cluster" ) {
#						name <- paste0("Cluster", name)
				}
				contrast <- strsplit(name, "-vs-")[[1]]
				dt <- stat2 %>% filter(.data[[group.by]] == i) %>% select(- !! group.by) %>%
					rename(!! contrast[[1]] := cells.1, !! contrast[[2]] := cells.2)
				WriteTable(dt, file = paste0("Stat.", name, ".xls"))
		}

		if ( do.plot ) {
				stat <- reshape2::melt(stat2, id.var = c("cluster", "contrast"), measure.vars = c("up", "down"), variable.name = "sig", value.name = "count" )
				plist <- list()
				for ( i in levels(stat[[group.by]]) ) {
						name <- i
						if ( group.by == "contrast" ) {
#								name <- paste0("Contrast.", name)
								group.by.inv <- "cluster"
						} else if ( group.by == "cluster" ) {
#								name <- paste0("Cluster", name)
								group.by.inv <- "contrast"
						}
						dt <- stat %>% filter(.data[[group.by]] == i) %>%
								droplevels() %>% 
								filter(count != 0)
						plist[[i]] <- ggplot(dt, aes_string(x = group.by.inv, y = "count", fill = "sig")) + 
								geom_bar(stat = "identity", position = "dodge") +
								geom_text(aes(label = count), hjust = 0.5, vjust = -0.5, position = position_dodge(0.9)) + 
								scale_fill_discrete(NULL, drop = FALSE) +
								scale_x_discrete(NULL, drop = FALSE) + 
								ylab("Number of Genes") + theme_light() +
								theme(axis.text.x = element_text(angle=45, hjust=1))
						plist[[i]] <- plist[[i]] + bar_theme_default()
#						x.text <- unique(p$data[[as.character(p$mapping$x)[2]]])
#						w <- grid::convertWidth(unit(1, "strwidth", x.text[which.max(nchar(x.text))]), "in", T) * length(x.text)
						w <- max(7, length(unique(stat[[group.by.inv]])) * 0.8)
						ggsave(plist[[i]], file = paste0("Stat.", name, ".pdf"), width = w, height = 7)
				}
		}
}


CalAvgByGrpInCls <- function(object, group.data, group_names,
				features = NULL, group.by = NULL, assay = NULL, slot = "data",
				is.expm1 = ifelse(slot=="data", TRUE, FALSE), diff_type = "all",
				outpref = "AllGene.avg_exp"
				){
		data <- GetAssayData(object, assay = assay, slot = slot)
		if ( ! is.null(features) ) data <- data[features, ]
		if ( is.expm1 ) data <- expm1(data)
		if ( ! is.null(group.by) ) Idents(object) <- group.by
		all_mean_exp <- NULL
		clusters <- if ( diff_type == "all" )  {
				c("bulk", levels(Idents(object)))
			} else if ( diff_type == "no_bulk" ) {
				levels(Idents(object))
			} else if ( diff_type == "only_bulk" ) {
				"bulk"
			} else {
				stop()
			}


		Gene_ID   <- ChangeOUTName(rownames(data), object@misc$fdata)
		Gene_name <- FindFeaturesName(object, rownames(data), "name")
		for ( i in clusters ) {
				cell.cluster <- if ( i == "bulk" ) !is.na(Idents(object)) else Idents(object) == i
				mean_exp <- do.call(cbind, 
								lapply(group_names, function(j){
										cell.group <- group.data[[j]]
										if ( sum(cell.cluster & cell.group) )
											Matrix::rowMeans(data[, cell.cluster & cell.group, drop = F])
										else
											rep(0, nrow(data))
								})
							)
				colnames(mean_exp) <- group_names
				if ( ! is.null(outpref) ) {						
						out_mean_exp <- cbind(Gene_ID = Gene_ID, Gene_name = Gene_name, mean_exp)
						WriteTable(out_mean_exp, file = paste(outpref, gsub("[/ ,()]", "_", i), "xls", sep = "."))
				}
				mean_exp <- data.frame(Cluster = i, mean_exp)
				all_mean_exp <- rbind(all_mean_exp, mean_exp)
		}
		WriteTable(clusters, file = "enrich.list", col.names = FALSE)

		invisible(all_mean_exp)
}


