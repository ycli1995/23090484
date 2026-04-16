
source("/state/partition1/NASDATA1/pipeline/Toolkit_ycli/Eula/v2.0/R/Eula.Seurat.R", chdir = TRUE)

df = ReadTable("diff_stat.xls")

df = df %>%
  filter(!grepl("\\-2", group)) %>%
  mutate(group = gsub("\\-1", "", group)) %>%
  mutate(group = factor(group, unique(group)))
WriteTable(df, "diff_stat_used.xls")

df = df %>%
  reshape2::melt("group")

colors = c("up" = "#d62926", "down" = "#2c77ac")
p = ggplot(df, aes(x = group, y = value)) +
  geom_bar(aes(fill = variable), stat = "identity", position = "dodge") +
  geom_text(aes(label = value, group = variable), vjust = -0.5, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Number of Metabolites", title = "Differential Metabolites Statistics") +
  ylim(0, 1.1 * max(df$value)) +
  bar_theme_default() +
  theme(axis.text  = element_text(color = "#000000", size = 12))
ggsave(p, filename = "diff_stat.pdf", width = 5.75, height = 4.75)

