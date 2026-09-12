suppressPackageStartupMessages({
  library(tidyverse)
})

# Paths are relative to the repository root; python run_all.py runs every module there.
data_dir <- file.path("intermediate", "clean")
out_dir <- file.path("intermediate", "analysis")
if (!dir.exists(data_dir)) stop("Run from the repository root after code/load_data.py has built intermediate/clean/", call. = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(0.10, "cm"),
      axis.title = element_text(size = axis_title_size, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.key = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

metric_labels <- c(
  laeq_24h_leq_db_a = "LAeq_24h",
  lday_leq_db_a = "Lday",
  lnight_leq_db_a = "Lnight",
  lpeak_leq_db_a = "Lpeak",
  linter_peak_leq_db_a = "Linter-peak",
  mp_leq_db_a = "Morning peak",
  ep_leq_db_a = "Evening peak",
  eve_op_leq_db_a = "Evening off-peak",
  laeq_sd_24h_db_a = "24h SD",
  laeq_peak_trough_diff_db = "Peak-trough",
  high_noise_share_gt_60db = "Share >60 dB",
  high_noise_share_gt_65db = "Share >65 dB"
)

noise <- read_csv(file.path(data_dir, "noise_metrics_station.csv"), show_col_types = FALSE)
master <- read_csv(file.path(data_dir, "station_master.csv"), show_col_types = FALSE) %>%
  select(station_key, station_name_cn, line_name, station_type_code, station_type_name_cn, day_type)

metrics <- intersect(names(metric_labels), names(noise))
df <- noise %>%
  filter(!is.na(laeq_24h_leq_db_a)) %>%
  select(station_key, all_of(metrics)) %>%
  drop_na(all_of(metrics)) %>%
  left_join(master, by = "station_key")

metric_mat <- df %>% select(all_of(metrics))
pearson_cor <- cor(metric_mat, use = "pairwise.complete.obs", method = "pearson")
spearman_cor <- cor(metric_mat, use = "pairwise.complete.obs", method = "spearman")

cor_long <- as.data.frame(as.table(spearman_cor)) %>%
  as_tibble() %>%
  rename(metric_x = Var1, metric_y = Var2, spearman_rho = Freq) %>%
  mutate(
    pearson_r = as.vector(pearson_cor),
    metric_x_label = recode(metric_x, !!!metric_labels),
    metric_y_label = recode(metric_y, !!!metric_labels)
  )
write_csv(cor_long, file.path(out_dir, "18_outcome_structure__metric_correlations.csv"))

pca <- prcomp(metric_mat, center = TRUE, scale. = TRUE)
pca_variance <- tibble(
  component = paste0("PC", seq_along(pca$sdev)),
  eigenvalue = pca$sdev^2,
  variance_explained = eigenvalue / sum(eigenvalue),
  cumulative_variance = cumsum(variance_explained)
)
write_csv(pca_variance, file.path(out_dir, "18_outcome_structure__pca_variance.csv"))

pca_loadings <- as_tibble(pca$rotation, rownames = "metric") %>%
  mutate(metric_label = recode(metric, !!!metric_labels)) %>%
  relocate(metric_label, .after = metric)
write_csv(pca_loadings, file.path(out_dir, "18_outcome_structure__pca_loadings.csv"))

station_scores <- as_tibble(pca$x[, 1:min(4, ncol(pca$x))], .name_repair = "minimal") %>%
  set_names(paste0("PC", seq_len(ncol(.)))) %>%
  bind_cols(df %>% select(station_key, station_name_cn, line_name, station_type_code, station_type_name_cn, day_type), .)
write_csv(station_scores, file.path(out_dir, "18_outcome_structure__station_scores.csv"))

summary_tbl <- tibble(
  item = c(
    "n_stations_complete",
    "n_metrics",
    "pc1_variance",
    "pc1_pc2_cumulative_variance",
    "mean_abs_spearman_offdiag",
    "min_spearman_offdiag",
    "max_spearman_offdiag"
  ),
  value = c(
    nrow(df),
    length(metrics),
    pca_variance$variance_explained[1],
    pca_variance$cumulative_variance[2],
    mean(abs(spearman_cor[upper.tri(spearman_cor)])),
    min(spearman_cor[upper.tri(spearman_cor)]),
    max(spearman_cor[upper.tri(spearman_cor)])
  )
)
write_csv(summary_tbl, file.path(out_dir, "18_outcome_structure__summary.csv"))

heatmap_df <- cor_long %>%
  mutate(
    metric_x_label = factor(metric_x_label, levels = rev(metric_labels[metrics])),
    metric_y_label = factor(metric_y_label, levels = metric_labels[metrics])
  )

p_heat <- ggplot(heatmap_df, aes(metric_y_label, metric_x_label, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-1, 1), name = "rho") +
  labs(x = NULL, y = NULL) +
  coord_equal() +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )
ggsave(file.path(out_dir, "18_outcome_structure__metric_correlation_heatmap.png"), p_heat, width = 178 / 25.4, height = 150 / 25.4, dpi = 600, bg = "white")

p_scree <- ggplot(pca_variance %>% slice(1:8), aes(component, variance_explained)) +
  geom_col(fill = "#0072B2", width = 0.7) +
  geom_line(aes(y = cumulative_variance, group = 1), colour = "#D55E00", linewidth = 0.5) +
  geom_point(aes(y = cumulative_variance), colour = "#D55E00", size = 1.4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Principal component", y = "Variance explained") +
  theme_pub()
ggsave(file.path(out_dir, "18_outcome_structure__pca_scree.png"), p_scree, width = 85 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

message("Completed outcome-structure PCA for ", nrow(df), " stations and ", length(metrics), " metrics.")
