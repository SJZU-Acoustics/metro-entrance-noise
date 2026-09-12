# =============================================================================
# Figure 2 -- Citywide temporal acoustic burden
#   (a) Hourly mean and median LAeq with inter-quartile ribbon; diurnal trough
#       and peak marked.
#   (b) Station-by-hour LAeq matrix, stations ordered by LAeq_24h.
#   (c) Distribution of station-level LAeq_24h and Lnight.
# =============================================================================
source(file.path("code", "display", "style.R"))

hc <- read_csv(file.path(EXPLORE_OUT, "03_descriptive__hourly_curve_stats.csv"),
               show_col_types = FALSE)
shp <- read_csv(file.path(EXPLORE_CLEAN, "station_hour_panel.csv"), show_col_types = FALSE)
nm  <- read_csv(file.path(EXPLORE_CLEAN, "noise_metrics_station.csv"), show_col_types = FALSE)

BS <- 8; AT <- 9  # 3-in-a-row typography

# ---- Panel a: diurnal curve --------------------------------------------------
trough <- hc %>% slice_min(mean, n = 1)
peak   <- hc %>% slice_max(mean, n = 1)
amp    <- round(peak$mean - trough$mean, 2)
write_csv(hc, file.path(DATA_DIR, "fig2_hourly_curve.csv"))

pa <- ggplot(hc, aes(hour)) +
  geom_ribbon(aes(ymin = p25, ymax = p75), fill = "#0072B2", alpha = 0.15) +
  geom_line(aes(y = median, linetype = "Median"), colour = "#0072B2", linewidth = 0.5) +
  geom_line(aes(y = mean,   linetype = "Mean"),   colour = "#0072B2", linewidth = 0.7) +
  geom_point(data = bind_rows(trough, peak), aes(y = mean), size = 1.5,
             colour = "#D55E00") +
  scale_linetype_manual(values = c(Mean = "solid", Median = "22"), name = NULL,
                        breaks = c("Mean", "Median")) +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24), expand = c(0.02, 0)) +
  scale_y_continuous(breaks = seq(52, 70, 3), limits = c(52, 70)) +
  labs(x = "Hour of day", y = expression(italic(L)[plain("Aeq")]~"(dB(A))")) +
  theme_pub(BS, AT) +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.06),
        legend.justification = c(1, 0), legend.key.width = unit(0.55, "cm"),
        legend.text = element_text(size = 7))

# ---- Panel b: station-hour matrix --------------------------------------------
ord <- nm %>% filter(!is.na(laeq_24h_leq_db_a)) %>% arrange(laeq_24h_leq_db_a) %>%
  mutate(rank = row_number()) %>% select(station_key, rank, laeq_24h_leq_db_a)

hm <- shp %>% filter(station_key %in% ord$station_key) %>%
  select(station_key, hour, laeq_hour_db_a) %>%
  left_join(ord, by = "station_key")
write_csv(hm, file.path(DATA_DIR, "fig2_station_hour_matrix.csv"))

pb <- ggplot(hm, aes(hour, rank, fill = laeq_hour_db_a)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B", name = expression(italic(L)[plain("Aeq")]~"(dB(A))"),
                       guide = guide_colourbar(barwidth = unit(0.30, "cm"),
                                               barheight = unit(2.4, "cm"))) +
  scale_x_continuous(breaks = seq(0, 24, 6), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(1, 109), expand = c(0, 0)) +
  labs(x = "Hour of day",
       y = expression("Station (ranked by "*italic(L)[plain("Aeq")]*" 24 h)")) +
  theme_pub(BS, AT) +
  theme(legend.position = "right", legend.title = element_text(size = 7, angle = 90),
        legend.title.position = "left",
        axis.line = element_blank(), axis.ticks = element_line(linewidth = 0.3))

# ---- Panel c: outcome distributions ------------------------------------------
dist <- nm %>% filter(!is.na(laeq_24h_leq_db_a)) %>%
  select(station_key, LAeq_24h = laeq_24h_leq_db_a, Lnight = lnight_leq_db_a) %>%
  pivot_longer(-station_key, names_to = "metric", values_to = "level")
write_csv(dist, file.path(DATA_DIR, "fig2_distributions.csv"))
means <- dist %>% group_by(metric) %>% summarise(m = mean(level), .groups = "drop")

# Outcome palette, kept consistent with Fig 5b (LAeq blue, Lnight purple).
DCOL <- c(LAeq_24h = "#0072B2", Lnight = "#CC79A7")
pc <- ggplot(dist, aes(level, fill = metric, colour = metric)) +
  geom_density(alpha = 0.30, linewidth = 0.5, adjust = 1.1) +
  geom_vline(data = means, aes(xintercept = m, colour = metric),
             linetype = "22", linewidth = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = DCOL, name = NULL,
                    labels = c(expression(italic(L)[plain("Aeq")]~"24 h"),
                               expression(italic(L)[plain("night")]))) +
  scale_colour_manual(values = DCOL, guide = "none") +
  scale_x_continuous(breaks = seq(45, 75, 10)) +
  labs(x = "Sound level (dB(A))", y = "Density") +
  theme_pub(BS, AT) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.text = element_text(hjust = 0))

# ---- Compose -----------------------------------------------------------------
fig2 <- pa + pb + pc +
  plot_layout(widths = c(1.15, 1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY))

save_fig(fig2, "fig2_temporal_burden.png", width_mode = "double_column", height_in = 2.7)
message("Fig 2 done. amplitude = ", amp, " dB")
