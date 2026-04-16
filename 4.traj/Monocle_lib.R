
CreateCDS <- function(
    obj, 
    parameter,
    assay = parameter$assay,
    slot = parameter$slot,
    run.disp = TRUE,
    feature_use = parameter$feature_use,
    feature_method = c("default", "Seurat", "monocle2"),
    nfeatures = 2000,
    mean_expression = 0.1,
    dispersion_empirical = 1,
    mnc_obj = parameter$mnc_obj,        
    ...
) {
  if (!is.null(mnc_obj) && file.exists(mnc_obj)) {
    return(Load(mnc_obj))
  }
  feature_method <- match.arg(feature_method)
  obj <- obj %>% 
    RenameSeuratColumns(parameter$rename_colnames) %>%
    RenameSeuratWrapper(parameter$rename) %>%
    SubsetObjectWrapper(parameter$subset) %>%
    CheckDefaultColumns(parameter$default_colnames)
  assay <- assay %||% DefaultAssay(obj)
  slot <- slot %||% "counts"
  data <- GetAssayData(obj[[assay]], slot)
  featuresID <- getFeaturesID(obj, rownames(data))
  data <- data[featuresID, , drop = FALSE]
  rownames(data) <- ChangeOUTName(rownames(data), obj@misc$fdata)
  # variable features
  order_features <- NULL
  if (length(feature_use) > 0) {
    order_features <- read.table(
      feature_use, 
      header = F, 
      stringsAsFactors = F
    )[, 1]
    feature_method <- ""
  }
  if (feature_method == "default") {
    order_features <- VariableFeatures(obj)
  }
  if (feature_method == "Seurat") {
    obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
    order_features <- VariableFeatures(obj)
  }
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }
  if (all(is.wholenumber(data[, head(1:ncol(data), 10000)]))) {
    expressionFamily <- negbinomial.size()
  } else if (any(data < 0)) {
    expressionFamily <- uninormal()
  } else {
    expressionFamily <- tobit()
  }
  
  gene_name <- getFeaturesName(obj, rownames(data))
  fd <- data.frame(gene_short_name = gene_name, row.names = rownames(data)) %>%
    new(Class = "AnnotatedDataFrame")
  pd <- obj[[]][colnames(data), , drop = F] %>% 
    new(Class = "AnnotatedDataFrame")
  valid_data <- data[, rownames(pd)]
  
  lowerDetectionLimit <- 0
  mnc_obj <- newCellDataSet(
    data, 
    phenoData = pd, 
    featureData = fd, 
    lowerDetectionLimit = lowerDetectionLimit, 
    expressionFamily = expressionFamily
  )
  mnc_obj <- estimateSizeFactors(mnc_obj)
  if (!run.disp) {
    return(mnc_obj)
  }
  mnc_obj <- estimateDispersions(mnc_obj)
  if (length(order_features) == 0) {
    # monocle2's ordering features
    negbin.names <- c('negbinomial', 'negbinomial.size')
    if (!any(familyname(mnc_obj@expressionFamily) %in% negbin.names)) {
      stop(familyname(mnc_obj@expressionFamily))
    }
    min.mean_expression <- mean_expression
    min.dispersion_empirical <- dispersion_empirical
    order_features <- dispersionTable(mnc_obj) %>%
      subset(
        mean_expression >= min.mean_expression &
          dispersion_empirical >= min.dispersion_empirical
      ) %>%
      pull("gene_id")
  }
  order_features <- getFeaturesID(obj, order_features)
  order_features <- ChangeOUTName(order_features, obj@misc$fdata)
  message("Ordering features (", length(order_features), "):")
  message(paste(head(order_features), collapse = ", "))
  mnc_obj <- setOrderingFilter(mnc_obj, ordering_genes = order_features)
  mnc_obj
}

select_ncenter <- function(mnc_obj, norm_method = c("log", "vstExprs", "none")) {
  if (ncol(mnc_obj) >= 100) {
    ncenter <- monocle:::cal_ncenter(ncol(mnc_obj))
  } else {
    ncenter <- ncol(mnc_obj) - 1
  }
  FM <- monocle:::normalize_expr_data(mnc_obj, norm_method, 1)
  xm <- Matrix::rowMeans(FM)
  xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
  FM <- FM[xsd > 0, , drop = FALSE]
  FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
  FM <- FM[!is.na(row.names(FM)), , drop = FALSE]
  FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), , drop = FALSE]
  X <- FM
  W <- DDRTree:::pca_projection_R(X %*% t(X), 2)
  Z <- t(W) %*% X
  select <- !sapply(
    X = seq(ncenter, 0), 
    FUN = function(i) t(Z)[seq(1, ncol(Z), length.out = i), ] %>% 
      duplicated() %>%
      any()
  )
  ncenter <- seq(ncenter, 0)[select][1]
  return(ncenter)
}

DoMonocle <- function( mnc_obj, ordering_genes = NULL, ... ) {
  #				max_components = 2, maxIter = 20, ncenter = NULL, tol = 0.001, sigma = 0.001, lambda = NULL, param.gamma = 10 ) {
  if (!is.null(ordering_genes)) {
    mnc_obj <- setOrderingFilter(mnc_obj, ordering_genes = ordering_genes)
  }
  if (!("negbinomial" == familyname(mnc_obj@expressionFamily)) || ("negbinomial.size" == familyname(mnc_obj@expressionFamily)) ) {
    norm_method <- "none"
  } else {
    norm_method <- c("log", "vstExprs", "none")
  }
  ncenter <- select_ncenter(mnc_obj, norm_method = norm_method)
  if (length(ncenter) == 0) {
    stop('unable to select proper ncenter, change ordering_genes may help.')
  }
  mnc_obj <- reduceDimension(mnc_obj, 
                             reduction_method = 'DDRTree', 
                             verbose = F, 
                             scaling = T, 
                             norm_method = norm_method, 
                             ncenter = ncenter, 
                             ...)
  mnc_obj <- orderCells(mnc_obj)
  return(mnc_obj)
}

RunMonocle <- function(
    mnc_obj, 
    parameter = list(), 
    feature_use = parameter$feature_use, 
    ...
) {
  genes.use <- if (!is.null(feature_use) && file.exists(feature_use[1])) {
    readLines(feature_use)
  }  else {
    NULL
  }
  mnc_obj <- DoMonocle(mnc_obj, ordering_genes = genes.use, ...)
  pData(mnc_obj) <- droplevels(pData(mnc_obj))
  return(mnc_obj)
}

RenameState <- function(
    mnc_obj, 
    parameter = list(), 
    rename_state = parameter$rename_state, 
    root_state = parameter$root_state, 
    reverse = parameter$reverse, 
    force = parameter$force
) {
  force <- force %||% FALSE
  old.state <- pData(mnc_obj) %>%
    top_n(-1, wt = Pseudotime) %>%
    pull("State")
  if (!is.null(root_state) ) {
    if (isTRUE(reverse)) {
      force <- TRUE  ## 暂时用来单独处理reverse的情况，后续版本更新再考虑重构逻辑
    }
    if (force || old.state != root_state) {
      if (isTRUE(reverse)) {
        root_state <- NULL
        message("------- doing reverse")
      }
      mnc_obj <- orderCells(mnc_obj, root_state = root_state, reverse = reverse)
    }
  }
  if (is.null(rename_state)) {
    return(mnc_obj)
  }
  if (length(rename_state) < length(levels(pData(mnc_obj)$State))) {
    stop("length of <rename_state> is less than State's levels.") 
  }
  levels(pData(mnc_obj)$State) <- rename_state
  pData(mnc_obj)$State <- factor(as.character(pData(mnc_obj)$State), levels = stringr::str_sort(levels(pData(mnc_obj)$State), numeric = T))
  print(str(pData(mnc_obj)))
  return(mnc_obj)
}

findRC <- function(n) {
  row <- ceiling(sqrt(n))
  col <- ceiling(n/row)
  return(c(row, col))
}

PlotMonocle <- function(
    mnc_obj,
    color_by,
    colors = NULL,
    show_branch_points = FALSE, 
    scale.size = 6,
    outpfx = NULL,
    ...
) {
  facet.scale <- 0.7
  p1 <- plot_cell_trajectory(
    mnc_obj, 
    color_by = color_by, 
    show_branch_points = show_branch_points
  )
  if (length(colors) > 0) {
    if (is.numeric(p1$data[[color_by]])) {
      p1 <- p1 + scale_color_gradientn(colors = colors)
    } else {
      p1 <- p1 + scale_color_manual(values = colors)
    }
  }
  if (!is.numeric(p1$data[[color_by]])) {
    ncol <- ceiling(length(unique(p1$data[[color_by]])) / 12)
    p1 <- p1 +
      guides(color = guide_legend(ncol = ncol)) +
      theme(legend.position = "right", legend.direction = "vertical")
    rc <- pData(mnc_obj)[[color_by]] %>% 
      levels() %>% 
      length() %>% 
      findRC()
    p2 <- p1 + facet_wrap(paste("~", color_by), nrow = rc[1])
    p <- list(p1, p2)
  } else {
    p <- list(p1)
  }
  if (length(outpfx) == 0) {
    return(p)
  }
  if (is.numeric(p1$data[[color_by]])) {
    lgd.size <- 0.8
  }	else {
    lgd.size <- max(1, max(nchar(pData(mnc_obj)[[color_by]] %>% as.character())) * 0.075) * ncol
  }
  ggsave(paste0(outpfx, ".pdf"), p[[1]], 
         width = scale.size + lgd.size, 
         height = scale.size, 
         limitsize = FALSE)
  if (length(p) > 1) {
    ggsave(paste0(outpfx, ".facet.pdf"), p[[2]], 
           width = scale.size * rc[2] * facet.scale + lgd.size, 
           height = scale.size * rc[1] * facet.scale, 
           limitsize = FALSE)
  }
  return(invisible(NULL))
}

select_sig_gene <- function(data, threshold = 1e-5, n = NULL, use.q = TRUE) {
  if (use.q) {
    tmp <- subset(data[order(data$qval),], qval < threshold )
  } else {
    tmp <- subset( data[order(data$pval),], pval < threshold )
  }
  if (is.null(n)) {
    return( rownames(tmp) )
  }
  return( head(rownames(tmp), n) )
}

FindBranchIndex <- function( cds ) {
  if (cds@dim_reduce_type == "DDRTree" ){
    reduced_dim_coords <- reducedDimK(cds)
  } else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  ica_space_df <- as.data.frame(Matrix::t(reduced_dim_coords)) %>%
    mutate(pseudo_point = rownames(.))
  mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
  branch_point_df <- ica_space_df %>%
    slice(match(mst_branch_nodes, pseudo_point)) %>%
    mutate(branch_point_idx = seq_len(n())) 
  BranchIndex <- branch_point_df$branch_point_idx %>%
    setNames(branch_point_df$pseudo_point)
  return(BranchIndex)
}

GetCellsInPath <- function(cds, branch_point = NULL, root_cell = NULL) {
  if (is.null(root_cell)) {
    root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  }
  if (cds@dim_reduce_type == "DDRTree") {
    pr_graph_cell_proj_mst <- minSpanningTree(cds)
  } else {
    pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
  }
  closest_vertex <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex
  vertex2cell <- by(
    rownames(closest_vertex), 
    paste0("Y_", closest_vertex), 
    function(x) as.character(x)
  )
  root_cell_Y <- paste0("Y_", closest_vertex[root_cell, 1])
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, branch_point, root_cell_Y)
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
  mst_no_branch_point <- pr_graph_cell_proj_mst - 
    V(pr_graph_cell_proj_mst)[branch_point]
  CellsInPath <- list()
  for (backbone_nei in neighbors(pr_graph_cell_proj_mst, branch_point)$name) {
    descendents <- bfs(
      mst_no_branch_point, 
      V(mst_no_branch_point)[backbone_nei], 
      unreachable = FALSE
    )
    descendents <- descendents$order[!is.na(descendents$order)]
    descendents <- V(mst_no_branch_point)[descendents]$name
    if (root_cell_Y %in% descendents == FALSE) {
      path_to_point <- c(rev(path_to_ancestor), branch_point, descendents ) %>%
        unique()
      if (cds@dim_reduce_type == "DDRTree") {
        cells_in_path <- unlist(vertex2cell[path_to_point], use.names = F)
      } else {
        cells_in_path <- path_to_point
      }
      cells_in_path <- intersect(cells_in_path, colnames(cds))
      CellsInPath[[backbone_nei]] <- cells_in_path
    }
  }
  return(CellsInPath)
}

.SummeriseCellsProperty <- function(
    cds, 
    cells = NULL, 
    col_name = 'State', 
    collapse = TRUE
) {
  if (is.null(cells)) {
    cells <- colnames(cds)
  }
  sub_p <- pData(cds)[cells, col_name, drop = F]
  p_name <- sub_p[[col_name]] %>% 
    unique() %>% 
    as.character()
  if (!collapse) {
    return(p_name)
  }
  name <- paste0("State ", paste(p_name, collapse = ","))
  return(name)
}

FindBranchName <- function(cds, branch_point = NULL, root_cell = NULL) {
  cells_in_path <- GetCellsInPath(cds, branch_point, root_cell)
  BranchName <- lapply(cells_in_path, 
                       .SummeriseCellsProperty, 
                       cds = cds, 
                       col_name = "State")
  return(BranchName)
}

ChangeBranchPointState <- function(cds) {
  if (cds@dim_reduce_type == "DDRTree") {
    pr_graph_cell_proj_mst <- minSpanningTree(cds)
  } else {
    pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
  }
  
  closest_vertex      <- as.data.frame(cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_closest_vertex)
  closest_vertex$name <- paste0( "Y_", closest_vertex$V1 )
  root_cell      <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root_vertex    <- closest_vertex[root_cell, 2]
  
  for ( branch_point in cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points ) {
    path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, branch_point, root_vertex)
    path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
    for ( i in 2:length(path_to_ancestor) ) {
      next_cells_to_root <- rownames(subset(closest_vertex, name == path_to_ancestor[i]))
      if ( length(next_cells_to_root) != 0 ) break
    }
    next_state_to_root <- unique(pData(cds)[next_cells_to_root, "State"])
    if ( length(next_state_to_root) == 1 ) {
      branch_cell <- rownames(subset(closest_vertex, name == branch_point))
      pData(cds)[branch_cell, "State"] <- next_state_to_root
    } else if ( length(next_state_to_root) == 0 ){
    } else {
      stop(next_state_to_root)
    }
  }
  
  return(cds)
}


############### __subroutine__ #################
################################################

RunMonocle <- function(mnc_obj, parameter = list(), feature_use = parameter$feature_use, ...){
  genes.use <- if ( ! is.null(feature_use) && file.exists(feature_use) ) readLines(feature_use) else NULL
  mnc_obj <- DoMonocle( mnc_obj, ordering_genes = genes.use, ... )
  pData(mnc_obj) <- droplevels(pData(mnc_obj))
  return(mnc_obj)
}



StatTrajectory <- function(
    mnc_obj, 
    name = "Trajectory.", 
    scale.size = 6, 
    show_branch_points = FALSE
){  
  Trajectory.data <- pData(mnc_obj) %>%
    tibble::rownames_to_column(var = "cells.name") %>%
    select(cells.name, Samples, Clusters, Groups, Pseudotime, State)
  WriteTable(Trajectory.data, file = paste0(name, "data.xls"))
}

Diff.State <- function(
    mnc_obj, 
    thres = 1e-7, 
    cores = 1, 
    use.q = TRUE, 
    diff_state_res = NULL,
    n_clusters = NULL
) {
  if (is.null(diff_state_res)) {
    if (file.exists("diff_state.rds")) {
      diff_state_res <- readRDS("diff_state.rds")
    } else {
      diff_state_res <- differentialGeneTest(
        mnc_obj, 
        fullModelFormulaStr = "~State", 
        cores = cores
      )
    }
  }
  saveRDS(diff_state_res, file = "diff_state.rds")
  
  thres <- as.numeric(thres)
  sig_gene_names <- select_sig_gene(diff_state_res, thres, use.q = use.q)
  
  if (length(sig_gene_names) == 0) {
    return(invisible(NULL))
  }
  p1 <- plot_genes_jitter(
    mnc_obj[head(sig_gene_names, 10), ], 
    grouping = "State", 
    color_by = "State", 
    ncol = 5
  )
  ggsave("Diff.state.pdf", p1, width = 12, height = 5, limitsize = FALSE)
  
  ## for online report
  nstates <- n_clusters %||% max(5, length(levels(pData(mnc_obj)$State)))
  nstates <- min(nstates, length(head(sig_gene_names, 1000)))
  p3 <- plot_pseudotime_heatmap(
    mnc_obj[head(sig_gene_names, 1000), ], 
    num_clusters = nstates, 
    cores = cores, 
    show_rownames = T, 
    return_heatmap = T
  )
  ggsave(p3[["ph_res"]], 
         file = "Diff.pseudotime_heatmap.pdf", 
         width = 7, 
         height = min(30, max(7, length(sig_gene_names) * 0.1)), 
         limitsize = FALSE) ## file name !
  WriteTable(sig_gene_names, 
             file = "Diff.state_heatmap.list.xls", 
             col.names = F)
  saveRDS(p3, file = "Diff.pseudotime_heatmap.Rds") 
  tmp <- by(
    rownames(pData(mnc_obj)), 
    pData(mnc_obj)$State, 
    function(x) Matrix::rowMeans(Biobase::exprs(mnc_obj)[sig_gene_names, x, drop = F])
  )
  tmp <- do.call(cbind.data.frame, tmp)
  colnames(tmp) <- paste0("State", colnames(tmp))

  diff_state_sig <- diff_state_res[sig_gene_names, ]
  diff_state_sig <- data.frame(
    GeneID = rownames(diff_state_sig), 
    tmp, 
    P_value = diff_state_sig$pval, 
    FDR = diff_state_sig$qval, 
    row.names = rownames(diff_state_sig) 
  )
  cluster <- cutree(p3[["ph_res"]]$tree_row, nstates)[p3[["ph_res"]]$tree_row$order]
  diff_state_sig <- diff_state_sig[intersect(rownames(diff_state_sig), names(cluster)), ]
  diff_state_sig$Gene_cluster <- cluster[rownames(diff_state_sig)]
  diff_state_sig$Current_gene_cluster <- diff_state_sig$Gene_cluster
  WriteTable(diff_state_sig, file = "Diff.state.xls")
#  WriteTable(diff_state_sig, file = "Diff.pseudotime_heatmap.xls")
  invisible(diff_state_res[sig_gene_names, ])
}

Diff.Pseudotime <- function(
    mnc_obj,
    thres = 1e-7, 
    cores = 1, 
    use.q = TRUE, 
    diff_Pseudotime_res = NULL,
    n_clusters = NULL
) {
  if (is.null(diff_Pseudotime_res)) {
    if (file.exists("diff_Pseudotime.rds")) {
      diff_Pseudotime_res <- readRDS("diff_Pseudotime.rds")
    } else {
      diff_Pseudotime_res <- differentialGeneTest(
        mnc_obj, 
        fullModelFormulaStr = "~sm.ns(Pseudotime)", 
        cores = cores
      )
      diff_Pseudotime_res <- diff_Pseudotime_res[rownames(mnc_obj), ]
      diff_Pseudotime_res <- diff_Pseudotime_res[!is.na(diff_Pseudotime_res[[1]]), ]
    }
  }
  saveRDS(diff_Pseudotime_res, file = "diff_Pseudotime.rds")
  
  thres <- as.numeric(thres)
  sig_gene_names <- select_sig_gene(diff_Pseudotime_res, 
                                    thres = thres, 
                                    use.q = use.q)
  
  if (length(sig_gene_names) == 0) {
    return(invisible(NULL))
  }
  p2 <- plot_genes_in_pseudotime(
    mnc_obj[head(sig_gene_names, 10), ], 
    ncol = 5, 
    color_by = "Samples",
    min_expr = 0.1
  )
  colors <- experimentData(mnc_obj)@other$colors[["Samples"]]
  if (length(colors) > 0) {
    p2 <- p2 + scale_color_manual(values = colors)
  }
  ggsave("Diff.genes_in_pseudotime.pdf", p2, 
         width = 12, 
         height = 5, 
         limitsize = FALSE)
    
  nstates <- n_clusters %||% max(5, length(levels(pData(mnc_obj)$State)))
  nstates <- min(nstates, length(head(sig_gene_names, 1000)))
  p3 <- plot_pseudotime_heatmap(
    mnc_obj[sig_gene_names, ], 
    num_clusters = nstates, 
    cores = cores, 
    show_rownames = T, 
    return_heatmap = T
  )
  ggsave(p3[["ph_res"]], file = "Diff.pseudotime_heatmap.pdf", 
         width = 7, 
         height = min(30, max(7, length(sig_gene_names) * 0.1)), 
         limitsize = FALSE)
  WriteTable(sig_gene_names, 
             file = "Diff.pseudotime_heatmap.list.xls", 
             col.names = F)
  saveRDS(p3, file = "Diff.pseudotime_heatmap.Rds")
 
  cluster <- cutree(p3[["ph_res"]]$tree_row, nstates)[p3[["ph_res"]]$tree_row$order]
  diff_Pseudotime_sig <- diff_Pseudotime_res[intersect(sig_gene_names, names(cluster)), ]
  diff_Pseudotime_sig <- data.frame(
    GeneID = rownames(diff_Pseudotime_sig), 
    Gene_cluster = cluster[rownames(diff_Pseudotime_sig)], 
    P_value = diff_Pseudotime_sig$pval, 
    FDR = diff_Pseudotime_sig$qval, 
    row.names = rownames(diff_Pseudotime_sig)
  )
  diff_Pseudotime_sig$Current_gene_cluster <- diff_Pseudotime_sig$Gene_cluster
  WriteTable(diff_Pseudotime_sig, file = "Diff.pseudotime.xls")
#  WriteTable(diff_Pseudotime_sig, file = "Diff.pseudotime_heatmap.xls")
  invisible(diff_Pseudotime_res[sig_gene_names, ])
}

Diff.Branch <- function(
    mnc_obj, 
    thres = 1e-7, 
    cores = 1, 
    use.q = TRUE, 
    BEAM_res = NULL
) {
  branch_point <- mnc_obj@auxOrderingData[[mnc_obj@dim_reduce_type]]$branch_points
  branch_point_index <- FindBranchIndex(mnc_obj)
  
  if (is.null(BEAM_res)) {
    if (file.exists("BEAM.rds")) {
      BEAM_res <- readRDS("BEAM.rds")
    } else {
      BEAM_res <- list()
    }
  }
  total_BEAM_sig <- list()
  for (i in seq(branch_point)) {
    cells_in_path <- GetCellsInPath(mnc_obj, branch_point[i])
    common_cell <- Reduce(intersect, cells_in_path)
    if (length(cells_in_path) < 2 || any(sapply(cells_in_path, function(i) length(setdiff(i, common_cell))) == 0)) {
      next
    }
    branch_name <- FindBranchName(ChangeBranchPointState(mnc_obj), branch_point[i])
    if (!i %in% names(BEAM_res)) {
      beam_res <- BEAM(
        mnc_obj, 
        branch_point = i, 
        cores = cores, 
        branch_labels = branch_name
      )
      BEAM_res[[as.character(i)]] <- beam_res
    } else {
      beam_res <- BEAM_res[[as.character(i)]]
    }
    saveRDS(beam_res, file = paste0("BEAM.", i, ".rds"))
    
    thres <- as.numeric(thres)
    sig_gene_names <- select_sig_gene(beam_res, thres = thres, use.q = use.q)
    
    if (length(sig_gene_names) > 0) {
      p4 <- plot_genes_branched_pseudotime(
        mnc_obj[head(sig_gene_names, 10), ], 
        branch_point = i, 
        color_by = "Samples", 
        ncol = 5, 
        branch_labels = branch_name, 
        min_expr = 0.1
      )
      colors <- experimentData(mnc_obj)@other$colors[["Samples"]]
      if (length(colors) > 0) {
        p4 <- p4 + scale_color_manual(values = colors)
      }
      ggsave(p4, file = paste0("Branch.", branch_point_index[branch_point[i]], ".genes_pseudotime.pdf"), width = 12, height = 5, limitsize = FALSE)
      
      num_clusters <- min(length(sig_gene_names), 5)
      p5 <- plot_genes_branched_heatmap(
        mnc_obj[sig_gene_names, ], 
        branch_point = i, 
        num_clusters = num_clusters, 
        cores = cores, 
        show_rownames = T, 
        return_heatmap=T, 
        branch_labels = branch_name
      )
      ggsave(
        p5[["ph_res"]], 
        file = paste0("Branch.", branch_point_index[branch_point[i]], ".genes_heatmap.pdf"), 
        width = 7,
        height = min(30, max(7, length(sig_gene_names) * 0.1)), 
        limitsize = FALSE
      )
      WriteTable(
        sig_gene_names, 
        file = paste0( "Branch.", branch_point_index[branch_point[i]], ".genes_heatmap.list.xls"),
        col.names = F
      )
      saveRDS(p5, file = paste0("Branch.", branch_point_index[branch_point[i]], ".genes_heatmap.Rds"))
      
      if (!is.null(p5$annotation_row)) {
        BEAM_sig <- beam_res[beam_res$gene_short_name %in% rownames(p5$annotation_row),,drop = FALSE]
        Gene_cluster <- p5$annotation_row[as.character(BEAM_sig$gene_short_name), "Cluster"]
      } else {
        BEAM_sig <- beam_res[sig_gene_names,,drop = FALSE]
        Gene_cluster <- '-'
      }
      branch_node <- branch_point_index[branch_point[i]]
      BEAM_sig <- data.frame(
        branch_node = branch_node,
        branch = paste(branch_name, collapse = " -vs- "),
        GeneID = rownames(BEAM_sig), 
        Gene_cluster = Gene_cluster, 
        P_value = BEAM_sig$pval, 
        FDR = BEAM_sig$qval, 
        row.names = rownames(BEAM_sig)
      )
      BEAM_sig$Current_gene_cluster <- BEAM_sig$Gene_cluster
      WriteTable(
        BEAM_sig, 
        file = paste0( "Branch.", branch_point_index[branch_point[i]], ".depended_gene.xls" )
      )
      total_BEAM_sig[[i]] <- BEAM_sig
    }
  }
  saveRDS(BEAM_res, file = "BEAM.rds")
  
  total_BEAM_sig <- do.call(rbind, total_BEAM_sig)
  WriteTable(total_BEAM_sig, file = "Branch.depended_gene.xls")
  
  invisible(total_BEAM_sig)
}


GetTrajectoryData <- function (
    cds, 
    x = 1, 
    y = 2, 
    theta = 0, 
    is.return = FALSE
) {
  #    requireNamespace("igraph")
  sample_state <- pData(cds)$State
  lib_info_with_pseudo <- pData(cds)
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- reducedDimS(cds)
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- reducedDimK(cds)
  } else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  ica_space_df <- Matrix::t(reduced_dim_coords) %>% 
    as.data.frame() %>% 
    select(prin_graph_dim_1 := !!x, prin_graph_dim_2 := !!y) %>% 
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  dp_mst <- minSpanningTree(cds)
  if (is.null(dp_mst)) {
    stop("You must first call orderCells() before using this function")
  }
  edge_df <- dp_mst %>% 
    igraph::as_data_frame() %>% 
    select(source = from, target = to) %>% 
    left_join(
      ica_space_df %>% 
        select(
          source = sample_name, 
          source_prin_graph_dim_1 = prin_graph_dim_1, 
          source_prin_graph_dim_2 = prin_graph_dim_2
        ), 
      by = "source"
    ) %>% 
    left_join(
      ica_space_df %>% 
        select(
          target = sample_name,
          target_prin_graph_dim_1 = prin_graph_dim_1,
          target_prin_graph_dim_2 = prin_graph_dim_2
        ), 
      by = "target"
    )
  data_df <- t(monocle::reducedDimS(cds)) %>% 
    as.data.frame() %>% 
    select(data_dim_1 := !!x, data_dim_2 := !!y) %>% 
    tibble::rownames_to_column("sample_name") %>% 
    mutate(sample_state) %>% 
    left_join(
      lib_info_with_pseudo %>% 
        tibble::rownames_to_column("sample_name"), 
      by = "sample_name"
    )
  if (theta!= 0) {
    return_rotation_mat <- function(theta) {
      theta <- theta/180 * pi
      matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
    }
    rot_mat <- return_rotation_mat(theta)
    cn1 <- c("data_dim_1", "data_dim_2")
    cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
    cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
    data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
    edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
    edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
  }
  
  branch_point_index <- FindBranchIndex(cds)
  edge_df <- edge_df %>% 
    mutate(
      source_branch_point_index = if_else(is.na(branch_point_index[source]), "NULL", as.character(branch_point_index[source])),
      target_branch_point_index = if_else(is.na(branch_point_index[target]), "NULL", as.character(branch_point_index[target]))
    )
  data_df <- data_df %>% 
    select(sample_name, 
           data_dim_1, 
           data_dim_2, 
           Samples, 
           Clusters, 
           Pseudotime, 
           State, 
           Size_Factor)
  if (is.return) {
    return(list(cells = data_df, bone = edge_df))
  }
  WriteTable(data_df, file = "Trajectory.cell.data.xls")
  WriteTable(edge_df, file = "Trajectory.bone.data.xls")
}

GetGenesInPseudotime <- function(
    cds_subset, 
    min_expr = NULL, 
    trend_formula = "~ sm.ns(Pseudotime, df=3)", 
    cores = 1, 
    is.return = FALSE, 
    outfile = NULL
) {
  outfile <- outfile %||% "Trajectory.genes.pseudotime.data.xls" 
  cds_exprs <- exprs(cds_subset)
  if (familyname(cds_subset@expressionFamily) %in% c("negbinomial", "negbinomial.size")) {
    cds_exprs_rel <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
  } else {
    cds_exprs_rel <- cds_exprs
  }
  cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  cds_exprs_rel <- reshape2::melt(round(as.matrix(cds_exprs_rel)))
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  colnames(cds_exprs_rel) <- c("f_id", "Cell", "expression_relative")
  cds_exprs <- cds_exprs %>% 
    merge(cds_exprs_rel) %>%
    merge(fData(cds_subset), by.x = "f_id", by.y = "row.names") %>%
    merge(pData(cds_subset), by.x = "Cell", by.y = "row.names") %>%
    mutate(f_id = as.character(f_id))
  
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
  model_expectation <- genSmoothCurves(
    cds_subset, 
    cores = cores, 
    trend_formula = trend_formula, 
    relative_expr = T, 
    new_data = new_data
  )
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- reshape2::melt(
    model_expectation, 
    varnames = c("f_id", "Cell"), 
    value.name = "expectation"
  )
  cds_exprs <- merge(cds_exprs, expectation)
  
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  
  cds_exprs <- cds_exprs %>% 
    select(Cell, 
           GeneID = f_id, 
           GeneName = gene_short_name, 
           expression, 
           expression_relative, 
           expectation, 
           Pseudotime, 
           State, 
           Samples, 
           Clusters)
  if (is.return) {
    return(cds_exprs)
  }
  WriteTable(cds_exprs, file = outfile)
}

GetGenesBranchedPseudotime <- function(
    cds, 
    branch_point = 1, 
    branch_labels = NULL, 
    method = "fitting", 
    min_expr = NULL, 
    trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch", 
    reducedModelFormulaStr = NULL, 
    is.return = FALSE, 
    cores = 1, 
    outfile = NULL, 
    ...
) {
  outfile <- outfile %||% "Trajectory.genes.branch.data.xls"
  branch_states <- NULL # don't know what it used for yet
  #    if (is.null(reducedModelFormulaStr) == FALSE) {
  #        pval_df <- branchTest(cds, branch_states = branch_states, branch_point = branch_point, fullModelFormulaStr = trend_formula, reducedModelFormulaStr = "~ sm.ns(Pseudotime, df=3)", ...)
  #        fData(cds)[, "pval"] <- pval_df[row.names(cds), "pval"]
  #    }
  if ("Branch" %in% all.vars(terms(as.formula(trend_formula)))) {
    cds_subset <- buildBranchCellDataSet(
      cds = cds, 
      branch_states = branch_states, 
      branch_point = branch_point, 
      branch_labels = branch_labels, 
      progenitor_method = "duplicate", 
      ...
    )
  } else {
    cds_subset <- cds
    pData(cds_subset)$Branch <- pData(cds_subset)$State
  }
  
  cds_exprs <- exprs(cds_subset)
  if (familyname(cds_subset@expressionFamily) %in% c("negbinomial", "negbinomial.size")) {
    cds_exprs_rel <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
  } else {
    cds_exprs_rel <- cds_exprs
  }
  cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  cds_exprs_rel <- reshape2::melt(round(as.matrix(cds_exprs_rel)))
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  colnames(cds_exprs_rel) <- c("f_id", "Cell", "expression_relative")
  cds_exprs <- cds_exprs %>% 
    merge(cds_exprs_rel) %>% 
    merge(fData(cds_subset), by.x = "f_id", by.y = "row.names") %>%
    merge(pData(cds_subset), by.x = "Cell", by.y = "row.names") %>%
    mutate(Branch = as.factor(Branch))
  
  new_data <- data.frame(
    Pseudotime = pData(cds_subset)$Pseudotime, 
    Branch = pData(cds_subset)$Branch
  )
  full_model_expectation <- genSmoothCurves(
    cds_subset, 
    cores = cores, 
    trend_formula = trend_formula, 
    relative_expr = T, 
    new_data = new_data
  )
  colnames(full_model_expectation) <- colnames(cds_subset)
  expectation <- reshape2::melt(
    full_model_expectation, 
    varnames = c("f_id", "Cell"), 
    value.name = "full_model_expectation"
  )
  cds_exprs <- merge(cds_exprs, expectation)
  #    if (!is.null(reducedModelFormulaStr)) {
  #        reduced_model_expectation <- genSmoothCurves(cds_subset, cores = cores, trend_formula = reducedModelFormulaStr, relative_expr = T, new_data = new_data)
  #        colnames(reduced_model_expectation) <- colnames(cds_subset)
  #        cds_exprs$reduced_model_expectation <- apply(cds_exprs, 1, function(x) reduced_model_expectation[x[2], x[1]])
  #    }
  if (method == "loess") {
    cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
  cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr
  #    if (!is.null(reducedModelFormulaStr)) {
  #        cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
  #        cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
  #    }
  cds_exprs <- cds_exprs %>%
    mutate(State = as.factor(State)) %>%
    mutate(Branch = as.factor(Branch)) %>% 
    select(Cell, 
           GeneID = f_id, 
           GeneName = gene_short_name, 
           expression, 
           expression_relative, 
           expectation = full_model_expectation, 
           Pseudotime, 
           State, 
           Samples, 
           Clusters, 
           original_cell_id, 
           Branch)
  
  if (is.return) {
    return(cds_exprs)
  }
  WriteTable(cds_exprs, file = outfile)
}

PlotPseudotimeHeatmap <- function(
    mnc_obj, 
    features, 
    parameter, 
    cores = 1, 
    outfile = NULL,
    diff_res = NULL,
    position = parameter$legend$position,
    title = parameter$labs$title,
    color = unlist(parameter$legend$color),
    cluster_rows = parameter$hclust$row,
    scale.range = parameter$data$scale.range,
    xlab = parameter$labs$xlab, 
    ylab = parameter$labs$ylab,
    show_rownames = parameter$axis$y.text$show, 
    show_colnames = parameter$axis$x.text$show,
    fontsize = parameter$font$size, 
    rowname_type = parameter$data$rowname_type,
    height = parameter$height,
    width = parameter$width,
    ...
) {
  position <- position %||% "right"
  title <- title %||% "Heatmap"
  cluster_rows <- cluster_rows %||% TRUE
  if (cluster_rows) diff_res = NULL
  scale.range <- scale.range %||% 3
  xlab <- xlab %||% "pseudotime"
  ylab <- ylab %||% "Gene"
  show_rownames <- show_rownames %||% TRUE 
  show_colnames <- show_colnames %||% FALSE
  fontsize <- fontsize %||% 10 
  rowname_type <- rowname_type %||% "id"
  
  is.legend <- if(position == "none") FALSE else TRUE
  title <- if (is.null(title)) NA else title
  num_clusters <- min(length(levels(pData(mnc_obj)$State)), length(features))
  
  hmcols <- monocle:::blue2green2red(100)
  if (length(unlist(color)) >= 3) {
    hmcols <- colorRampPalette(color)(100)
  }
  ph <- plot_pseudotime_heatmap(
    mnc_obj[features, ], 
    num_clusters = num_clusters, 
    cores = cores, 
    show_rownames = T, 
    return_heatmap = T,
    use_gene_short_name = FALSE,
    cluster_rows = cluster_rows, 
    hmcols = hmcols, 
    scale_max = scale.range, 
    scale_min = scale.range * -1
  )
  height <- height %||% 
    grid::convertHeight(sum(ph$ph_res$gtable$heights), "inches", valueOnly = T)
  width  <- width %||% 
    grid::convertWidth(sum(ph$ph_res$gtable$widths), "inches", valueOnly = T)
  while (!is.null(dev.list())) { 
    dev.off()
  }
  
  if (!is.null(outfile)) {
    pdf(file = outfile, height = height, width = width)
  }
  if (!is.null(xlab) || !is.null(ylab)) {
    require(grid)
    setHook(
      "grid.newpage", 
      function() {
        pushViewport(viewport(
          x = 1, 
          y = 1, 
          width = 0.9, 
          height = 0.9, 
          name = "vp", 
          just = c("right","top")
        ))
      }, 
      action = "prepend"
    )
  }
  annotation_row <- ph$annotation_row
  mat <- ph$heatmap_matrix[, ]
  if (!is.null(diff_res)) {
    rownames(diff_res) <- diff_res$GeneID
    features <- intersect(features, diff_res$GeneID)
    annotation_row <- diff_res[features, , drop = FALSE] %>%
      mutate(Gene_cluster = diff_res[features, "Gene_cluster"]) %>%
      mutate(Gene_cluster = factor(Gene_cluster, stringr::str_sort(unique(Gene_cluster), numeric = TRUE))) %>%
      select(Gene_cluster) %>%
      arrange(Gene_cluster) %>%
      as.data.frame()
    num_clusters <- NA
    cluster_rows <- FALSE
  }
  if (!is.null(annotation_row)) {
    colnames(annotation_row) <- "Gene_cluster"
    mat <- mat[rownames(annotation_row), , drop = FALSE]
  }
  old.rownames <- rownames(mat)
  if ( rowname_type == "name" ) {
    name <- fData(mnc_obj)[rownames(mat), "gene_short_name"]
    rownames(mat) <- name
    if (!is.null(annotation_row)) rownames(annotation_row) <- name
  } 
  if (rowname_type == "id_name") {
    name <- fData(mnc_obj)[rownames(mat), "gene_short_name"]
    name <- paste0(rownames(mat), "(", name, ")")
    rownames(mat) <- name
    if (!is.null(annotation_row)) rownames(annotation_row) <- name
  }
  pheatmap::pheatmap(
    mat,
    useRaster = TRUE, 
    border_color = NA, 
    color = ph$hmcols, 
    breaks = ph$bks,
    clustering_method = ph$clustering_method,
    cutree_rows = num_clusters, 
    cluster_rows = cluster_rows, 
    cluster_cols = FALSE,
    clustering_distance_rows = ph$row_dist,
    treeheight_row = 20, 
    legend = is.legend,
    #annotation_row = ph$annotation_row, 
    annotation_row = annotation_row,
    annotation_col = ph$annotation_col,
    annotation_names_row = FALSE, 
    annotation_names_col = FALSE,
    show_rownames = show_rownames, 
    show_colnames = show_colnames,
    main = title,
    fontsize = fontsize,
    labels_row = NULL,
    ...
  )
  if (!is.null(xlab) || ! is.null(ylab)) {
    setHook("grid.newpage", NULL, "replace")
    grid.text(xlab, y = -0.05, gp = gpar(fontsize = fontsize))
    grid.text(ylab, x = -0.05, rot = 90, gp = gpar(fontsize = fontsize))
  }
  if (!is.null(outfile)) {
    dev.off()
  }
  return(list(ph = ph, annotation_row = annotation_row, old.rownames = old.rownames))
}

PlotBranchHeatmap <- function(
    mnc_obj, 
    features, 
    parameter, 
    cores = 1, 
    outfile = NULL, 
    branch_point = NULL,
    diff_res = NULL,
    position = parameter$legend$position,
    title = parameter$labs$title,
    color = unlist(parameter$legend$color),
    scale.range = parameter$data$scale.range,
    cluster_rows = parameter$hclust$row,
    xlab = parameter$labs$xlab, 
    ylab = parameter$labs$ylab,
    show_rownames = parameter$axis$y.text$show, 
    show_colnames = parameter$axis$x.text$show,
    fontsize = parameter$font$size, 
    rowname_type = parameter$data$rowname_type,
    height = parameter$height,
    width = parameter$width,
    ...
) {
  position <- position %||% "right"
  title <- title %||% "Heatmap"
  cluster_rows <- cluster_rows %||% TRUE
  if (cluster_rows) diff_res = NULL
  scale.range <- scale.range %||% 3
  xlab <- xlab %||% "pseudotime"
  ylab <- ylab %||% "Gene"
  show_rownames <- show_rownames %||% TRUE 
  show_colnames <- show_colnames %||% FALSE
  fontsize <- fontsize %||% 10 
  rowname_type <- rowname_type %||% "id"
  
  is.legend <- if(position == "none") FALSE else TRUE
  title <- if (is.null(title)) NA else title
  num_clusters <- min(length(levels(pData(mnc_obj)$State)), length(features))
  
  branch_point_index <- FindBranchIndex(mnc_obj)
  branch_point_name <- names(branch_point_index)[branch_point_index == branch_point]
  branch_point <- branch_point_index[branch_point_index == branch_point]
  branch_name <- FindBranchName(ChangeBranchPointState(mnc_obj), branch_point_name)
  
  hmcols <- monocle:::blue2green2red(100)
  if (length(unlist(color)) >= 3) {
    hmcols <- colorRampPalette(color)(100)
  }
  ph <- plot_genes_branched_heatmap(
    mnc_obj[features, ], 
    branch_point = branch_point,
    num_clusters = num_clusters, 
    cores = cores,
    show_rownames = T, 
    return_heatmap = T,
    branch_labels = branch_name, 
    cluster_rows = cluster_rows, 
    hmcols = hmcols,
    use_gene_short_name = FALSE,
    scale_max = scale.range,
    scale_min = scale.range * -1
  )
  height <- height %||% 
    grid::convertHeight(sum(ph$ph_res$gtable$heights), "inches", valueOnly = T)
  width  <- width %||% 
    grid::convertWidth(sum(ph$ph_res$gtable$widths), "inches", valueOnly = T)
  while (!is.null(dev.list())) { 
    dev.off()
  }
  
  if (!is.null(outfile)) {
    pdf(file = outfile, height = height, width = width)
  }
  if (!is.null(xlab) || !is.null(ylab)) {
    require(grid)
    setHook(
      "grid.newpage", 
      function() {
        pushViewport(viewport(
          x = 1, 
          y = 1, 
          width = 0.9, 
          height = 0.9, 
          name = "vp", 
          just = c("right","top")
        ))
      }, 
      action = "prepend"
    )
  }
  annotation_row <- ph$annotation_row
  mat <- ph$heatmap_matrix[, ]
  if (!is.null(diff_res)) {
    rownames(diff_res) <- diff_res$GeneID
    features <- intersect(features, diff_res$GeneID)
    diff_res <- diff_res[diff_res$branch_node == branch_point, , drop = FALSE]
    annotation_row <- diff_res[features, , drop = FALSE] %>%
      mutate(Gene_cluster = diff_res[features, "Gene_cluster"]) %>%
      mutate(Gene_cluster = factor(Gene_cluster, stringr::str_sort(unique(Gene_cluster), numeric = TRUE))) %>%
      select(Gene_cluster) %>%
      arrange(Gene_cluster) %>%
      as.data.frame()
    num_clusters <- NA
    cluster_rows <- FALSE
  }
  if (!is.null(annotation_row)) {
    colnames(annotation_row) <- "Gene_cluster"
    mat <- mat[rownames(annotation_row), , drop = FALSE]
  }
  old.rownames <- rownames(mat)
  if ( rowname_type == "name" ) {
    name <- fData(mnc_obj)[rownames(mat), "gene_short_name"]
    rownames(mat) <- name
    if (!is.null(annotation_row)) rownames(annotation_row) <- name
  } 
  if (rowname_type == "id_name") {
    name <- fData(mnc_obj)[rownames(mat), "gene_short_name"]
    name <- paste0(rownames(mat), "(", name, ")")
    rownames(mat) <- name
    if (!is.null(annotation_row)) rownames(annotation_row) <- name
  }
  
  pheatmap::pheatmap(
    mat,
    useRaster = T, 
    border_color = NA, 
    color = ph$hmcols, 
    breaks = ph$bks,
    clustering_method = "ward.D2",
    cutree_rows = num_clusters,
    cluster_rows = cluster_rows, 
    cluster_cols = FALSE,
    clustering_distance_rows = ph$row_dist,
    treeheight_row = 20, 
    legend = is.legend,
    annotation_row = annotation_row, 
    annotation_col = ph$annotation_col,
    annotation_colors = ph$annotation_colors,
    annotation_names_row = FALSE, 
    annotation_names_col = FALSE,
    show_rownames = show_rownames, 
    show_colnames = show_colnames,
    main = title,
    fontsize = fontsize,
    labels_row = NULL,
    gaps_col = ph$col_gap_ind,
    #				filename = outfile, width = 7, height = NA,
    ...
  )
  if (!is.null(xlab) || ! is.null(ylab)) {
    setHook("grid.newpage", NULL, "replace")
    grid.text(xlab, y = -0.05, gp = gpar(fontsize = fontsize))
    grid.text(ylab, x = -0.05, rot = 90, gp = gpar(fontsize = fontsize))
  }
  if (!is.null(outfile)) {
    dev.off()
  }
  return(list(ph = ph, annotation_row = annotation_row, old.rownames = old.rownames))
}

CDS_avg <- function(cds, group.by = "State", outfile = "AllGene.avg_exp.xls") {
  tmp <- by(
    rownames(pData(cds)), 
    pData(cds)[[group.by]], 
    function(x) Matrix::rowMeans(Biobase::exprs(cds)[, x, drop = F]) 
  )
  tmp <- do.call(cbind.data.frame, tmp)
  colnames(tmp) <- paste0(group.by, colnames(tmp))
  tmp <- data.frame(
    GeneID = rownames(tmp), 
    GeneName = fData(cds)$gene_short_name, 
    tmp
  )
  WriteTable(tmp, file = outfile )
}




