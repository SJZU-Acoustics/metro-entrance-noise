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
      legend.title = element_blank(),
      legend.key = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

cyclic_lag <- function(target_hour, reference_hour) {
  ((target_hour - reference_hour + 12) %% 24) - 12
}

safe_cor <- function(x, y, method) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 6 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

panel <- read_csv(file.path(data_dir, "station_hour_panel.csv"), show_col_types = FALSE) %>%
  mutate(
    laeq_hour_db_a = as.numeric(laeq_hour_db_a),
    pop_heat_50m = as.numeric(pop_heat_50m),
    pop_heat_20m = as.numeric(pop_heat_20m)
  )

valid_panel <- panel %>%
  filter(!is.na(laeq_hour_db_a), !is.na(pop_heat_50m))

station_sync <- valid_panel %>%
  group_by(station_key, station_name_cn, line_name, station_type_code, station_type_name_cn, day_type) %>%
  summarise(
    n_hours = n(),
    spearman_noise_heat_50m = safe_cor(laeq_hour_db_a, pop_heat_50m, "spearman"),
    pearson_noise_heat_50m = safe_cor(laeq_hour_db_a, pop_heat_50m, "pearson"),
    noise_peak_hour = hour[which.max(laeq_hour_db_a)][1],
    heat_peak_hour_50m = hour[which.max(pop_heat_50m)][1],
    noise_trough_hour = hour[which.min(laeq_hour_db_a)][1],
    heat_trough_hour_50m = hour[which.min(pop_heat_50m)][1],
    high_high_hours_q75 = sum(
      laeq_hour_db_a >= quantile(laeq_hour_db_a, 0.75, na.rm = TRUE) &
        pop_heat_50m >= quantile(pop_heat_50m, 0.75, na.rm = TRUE),
      na.rm = TRUE
    ),
    high_high_share_q75 = high_high_hours_q75 / n_hours,
    high_noise_high_heat_hours = sum(
      laeq_hour_db_a >= 65 &
        pop_heat_50m >= median(pop_heat_50m, na.rm = TRUE),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    peak_lag_heat_minus_noise_hours = cyclic_lag(heat_peak_hour_50m, noise_peak_hour),
    trough_lag_heat_minus_noise_hours = cyclic_lag(heat_trough_hour_50m, noise_trough_hour),
    synchrony_type = case_when(
      is.na(spearman_noise_heat_50m) ~ "insufficient",
      spearman_noise_heat_50m >= 0.50 ~ "strong positive",
      spearman_noise_heat_50m >= 0.20 ~ "moderate positive",
      spearman_noise_heat_50m > -0.20 ~ "weak or mixed",
      TRUE ~ "inverse"
    )
  )

write_csv(station_sync, file.path(out_dir, "19_temporal_synchrony__station_metrics.csv"))

sync_summary <- station_sync %>%
  summarise(
    n_stations = n(),
    median_spearman = median(spearman_noise_heat_50m, na.rm = TRUE),
    mean_spearman = mean(spearman_noise_heat_50m, na.rm = TRUE),
    p25_spearman = quantile(spearman_noise_heat_50m, 0.25, na.rm = TRUE),
    p75_spearman = quantile(spearman_noise_heat_50m, 0.75, na.rm = TRUE),
    n_positive_0_2 = sum(spearman_noise_heat_50m >= 0.20, na.rm = TRUE),
    n_positive_0_5 = sum(spearman_noise_heat_50m >= 0.50, na.rm = TRUE),
    n_inverse = sum(spearman_noise_heat_50m <= -0.20, na.rm = TRUE),
    median_abs_peak_lag = median(abs(peak_lag_heat_minus_noise_hours), na.rm = TRUE),
    mean_high_high_share_q75 = mean(high_high_share_q75, na.rm = TRUE)
  )
write_csv(sync_summary, file.path(out_dir, "19_temporal_synchrony__summary.csv"))

type_counts <- station_sync %>%
  count(synchrony_type, name = "n_stations") %>%
  mutate(share = n_stations / sum(n_stations))
write_csv(type_counts, file.path(out_dir, "19_temporal_synchrony__type_counts.csv"))

high_high_hours <- valid_panel %>%
  group_by(station_key) %>%
  mutate(
    high_high_q75 = laeq_hour_db_a >= quantile(laeq_hour_db_a, 0.75, na.rm = TRUE) &
      pop_heat_50m >= quantile(pop_heat_50m, 0.75, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(high_high_q75, !is.na(primary_event_cn))

event_comp <- high_high_hours %>%
  count(primary_event_cn, name = "n_hours") %>%
  arrange(desc(n_hours)) %>%
  mutate(share = n_hours / sum(n_hours))
write_csv(event_comp, file.path(out_dir, "19_temporal_synchrony__high_high_event_composition.csv"))

p_rho <- ggplot(station_sync, aes(spearman_noise_heat_50m)) +
  geom_histogram(binwidth = 0.1, boundary = 0, fill = "#0072B2", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = "Within-station Spearman rho", y = "Stations (n)") +
  theme_pub()
ggsave(file.path(out_dir, "19_temporal_synchrony__rho_distribution.png"), p_rho, width = 85 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

peak_df <- station_sync %>%
  count(noise_peak_hour, heat_peak_hour_50m, name = "n_stations")
p_peak <- ggplot(peak_df, aes(noise_peak_hour, heat_peak_hour_50m, size = n_stations)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, linetype = "dashed", colour = "grey40") +
  geom_point(colour = "#D55E00", alpha = 0.75) +
  scale_x_continuous(breaks = seq(0, 23, 4), limits = c(-0.5, 23.5)) +
  scale_y_continuous(breaks = seq(0, 23, 4), limits = c(-0.5, 23.5)) +
  scale_size_area(max_size = 5) +
  labs(x = "Noise peak hour", y = "Heat peak hour") +
  theme_pub() +
  theme(legend.position = "right")
ggsave(file.path(out_dir, "19_temporal_synchrony__peak_hour_alignment.png"), p_peak, width = 85 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

message("Completed temporal synchrony analysis for ", nrow(station_sync), " stations.")
