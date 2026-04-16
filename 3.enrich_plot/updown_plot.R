
library(ggplot2)
library(data.table)
library(dplyr)

diff = c("up", "down")

#go = list()
ko = list()
for (i in diff) {
  ko[[i]] = fread(paste0(i, ".bar_Gradient.xls"), data.table = F) %>%
    mutate(sig = i, type = "KEGG")
#  go[[i]] = fread(paste0(i, ".GO.bar_Gradient.xls"), data.table = F) %>%
#    mutate(sig = i, type = class)
}
#go = Reduce(rbind, go)
ko = Reduce(rbind, ko)

df = rbind(ko) %>%
  mutate(log10P = -log10(Pvalue)) %>%
  mutate(log10P = ifelse(sig == "up", log10P, -log10P)) %>%
  mutate(sig = factor(sig, diff))

max.x = max(abs(df$log10P))
cols = c("up" = "#ef6775", "down" = "#4077ab")
for (i in unique(df$type)) {
  i2 = gsub("\\s|\\/", "_", i)
  df2 = df %>% filter(type == i)

  for (top.n in c(5, 10)) {
    df3 = df2 %>%
      filter(Pvalue < 0.05) %>%
      group_by(sig) %>%
      slice_max(abs(log10P), n = top.n, with_ties = F) %>%
      ungroup() %>%
      mutate(Descrption = factor(Descrption, rev(unique(Descrption)))) %>%
      mutate(Descrption2 = paste(Descrption, sig)) %>%
      mutate(Descrption2 = factor(Descrption2, rev(unique(Descrption2))))
    p = ggplot(df3, aes(x = log10P, y = Descrption2)) +
      geom_bar(aes(fill = sig), stat = "identity") +
      scale_fill_manual(values = cols) + 
      geom_text(data = df3 %>% filter(sig == "up"), aes(label = Descrption, x = -0.02), hjust = 1) +
      geom_text(data = df3 %>% filter(sig == "down"), aes(label = Descrption, x = 0.02), hjust = 0) +
      xlim(-max.x, max.x) +
      labs(x = "-log10(Pvalue)", y = "", fill = "") +
      theme_classic() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    h = nrow(df3) * 0.275 + 1.6
    w = max(7, max(nchar(as.character(df3$Descrption))) * 0.08 * 2 + 0.75)
    print(c(h, w))
    ggsave(p, filename = paste0(i2, ".UpDown.", top.n, ".barplot.pdf"), width = w, height = h)
    system(paste("convert -density 300", paste0(i2, ".UpDown.", top.n, ".barplot.pdf"), paste0(i2, ".UpDown.", top.n, ".barplot.png")))
  }
}


