# SI Fig S1 -- Acoustic-outcome structure: PCA scree + metric correlation heatmap.
source(file.path("code", "display", "style.R"))

var <- read_csv(file.path(EXPLORE_OUT, "18_outcome_structure__pca_variance.csv"),
                show_col_types = FALSE)
mc  <- read_csv(file.path(EXPLORE_OUT, "18_outcome_structure__metric_correlations.csv"),
                show_col_types = FALSE)
BS <- 8; AT <- 9

lvl <- c("LAeq_24h","Lday","Lnight","Lpeak","Linter-peak","Morning peak",
         "Evening peak","Evening off-peak","24h SD","Peak-trough",
         "Share >60 dB","Share >65 dB")

# Panel a: scree (top 8 components)
va <- var %>% slice(1:8) %>%
  mutate(component = factor(component, levels = component))
pa <- ggplot(va, aes(component, variance_explained)) +
  geom_col(fill = "#0072B2", width = 0.7) +
  geom_line(aes(y = cumulative_variance, group = 1), colour = "grey30", linewidth = 0.5) +
  geom_point(aes(y = cumulative_variance), colour = "grey30", size = 1.1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Principal component", y = "Variance explained") +
  theme_pub(BS, AT) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.8))

# Panel b: 12x12 Spearman correlation heatmap
hm <- mc %>%
  mutate(x = factor(metric_x_label, levels = lvl),
         y = factor(metric_y_label, levels = rev(lvl))) %>%
  filter(!is.na(x), !is.na(y))
pb <- ggplot(hm, aes(x, y, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = DIVERGE_LOW, mid = DIVERGE_MID, high = DIVERGE_HIGH,
                       midpoint = 0, limits = c(-1, 1), name = "Spearman ρ",
                       guide = guide_colourbar(barwidth = unit(0.30, "cm"),
                                               barheight = unit(2.6, "cm"))) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_pub(BS, AT) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.2),
        axis.text.y = element_text(size = 6.2),
        axis.line = element_blank(), axis.ticks = element_line(linewidth = 0.3),
        legend.title = element_text(size = 7, angle = 90),
        legend.title.position = "left")

fig <- pa + pb + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY))
save_fig(fig, "si_outcome_structure.png", width_mode = "double_column", height_in = 3.3, dir = SIFIG_DIR)
message("SI S1 done.")
