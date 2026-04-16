
library(data.table)
library(dplyr)
library(ggplot2)

MeteorPlot <- function(
    df, 
    P.column = "Pvalue",
    desc.column = "Descrption",
    ratio.column = "ratio",
    counts.column = "num",
    use.counts = FALSE,
    facet.by = "class",
    colors = RColorBrewer::brewer.pal(4, "Set1"),
    top.n = 5,
    ...
) {
  df2 <- df %>%
    mutate(logFDR = -log10(!!sym(P.column))) %>%
    group_by(across(all_of(facet.by))) %>%
    arrange(desc(logFDR), .by_group = TRUE) %>%
    slice_max(n = top.n, order_by = logFDR, with_ties = FALSE) %>%
    mutate(Description = as.character(!!sym(desc.column))) %>%
    mutate(Description = factor(Description, rev(unique(Description))))
  
  max.logFDR = max(df2$logFDR)
  if (use.counts) {
    max.count = max(df2[[counts.column]])
  } else {
    max.count <- max(df2[[ratio.column]] * 100)
  }
  step.count = round(max.count / 5, digits = 2)
  break.count = seq(0, max.count, step.count)
  fold = max.count / max.logFDR
  
  print(max.count)
  print(break.count)
  print(fold)
  
  type.colors <- colors[seq_along(unique(df2[[facet.by]]))]
  names(type.colors) <- unique(df2[[facet.by]])
  
  p1 = ggplot(df2) +
    ggforce::geom_link(
      aes(
        x = 0, y = Description, 
        xend = logFDR, yend = Description, 
        color = !!sym(facet.by), 
        linewidth = after_stat(index)
      ),
      n = 500,
      show.legend = F
    ) +
    geom_point(aes(x = logFDR, y = Description), 
               color = "black", fill = "white", 
               size = 6, shape = 21)
  if (use.counts) {
    p1 = p1 + 
      geom_line(aes(x = !!sym(counts.column) / fold, y = Description, group = 1), 
                orientation = "y", linewidth = 1, color = "#FFCC00") +
      scale_x_continuous(sec.axis = sec_axis(~. * fold, name = "Gene number")) 
  } else {
    p1 = p1 + 
      geom_line(aes(x = !!sym(ratio.column) * 100 / fold, y = Descrption, group = 1), 
                orientation = "y", linewidth = 1, color = "#FFCC00") +
      scale_x_continuous(sec.axis = sec_axis(~. * fold, name = "Percent of geneRatio (%)")) 
  }
  p1 = p1 +
    facet_wrap(paste("~", facet.by), ncol = 1, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text = element_text(color = "black", size = 12), 
          plot.title = element_text(hjust = 0.5)) +
    labs(x = sprintf("-log10(%s)", P.column), y = "", title = paste("Enrichment")) +
    scale_color_manual(values = type.colors)
  
  return(p1)
}


EnrichBarplot <- function(
    df,
    P.column = "Pvalue",
    desc.column = "Descrption",
    ratio.column = "ratio",
    counts.column = "num",
    use.counts = FALSE,
    colors = RColorBrewer::brewer.pal(9, "YlOrRd")[1:8],
    top.n = 20,
    ...
) {
  df2 <- df %>%
    mutate(logFDR = -log10(!!sym(P.column))) %>%
    arrange(desc(logFDR)) %>%
    slice_max(n = top.n, order_by = logFDR, with_ties = FALSE) %>%
    mutate(Description = as.character(!!sym(desc.column))) %>%
    mutate(Description = factor(Description, rev(unique(Description))))
  
  if (use.counts) {
    y.column <- counts.column
    ylab <- "Gene number"
  } else {
    y.column <- ratio.column
    ylab <- "Gene ratio"
  }
  
  FDR.colors <- colors
  
  p1 = ggplot(df2, aes(y = !!sym(y.column), x = Description)) +
    geom_bar(aes(fill = logFDR), stat = "identity") +
    scale_fill_gradientn(colors = FDR.colors) +
    geom_text(aes(y = 0, label = Description), hjust = 0) +
    coord_flip() +
    labs(y = ylab, x = "Description", fill = sprintf("-log10(%s)", P.column), title = paste("Enrichment")) +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5))
  
  return(p1)
}

