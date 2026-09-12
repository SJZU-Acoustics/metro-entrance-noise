# =============================================================================
# Figure 3 -- Activity-synchronised noise dynamics within stations
#   (a) Distribution of within-station heat-noise Spearman rho.
#   (b) Alignment of the heat-activity peak hour and the noise peak hour.
#   (c) Sound-source composition of high-heat / high-noise hours.
# =============================================================================
source(file.path("code", "display", "style.R"))

sm  <- read_csv(file.path(EXPLORE_OUT, "19_temporal_synchrony__station_metrics.csv"),
                show_col_types = FALSE)
sumr <- read_csv(file.path(EXPLORE_OUT, "19_temporal_synchrony__summary.csv"),
                 show_col_types = FALSE)
comp <- read_csv(file.path(EXPLORE_OUT, "19_temporal_synchrony__high_high_event_composition.csv"),
                 show_col_types = FALSE)

BS <- 8; AT <- 9
SYN_COLS <- c("strong positive" = "#009E73", "moderate positive" = "#E69F00",
              "inverse" = "#D55E00")
SYN_LAB  <- c("strong positive" = "Strong (ρ ≥ 0.5)",
              "moderate positive" = "Moderate",
              "inverse" = "Inverse")
# Panel (c) source-group palette: a neutral, low-saturation blue-grey ramp that
# does not collide with the synchrony-class greens/oranges in panels (a)/(b)
# (SOURCE_COLS reuses green + vermillion, so it is reserved for Figs 4-5 only).
COMP_COLS <- c(Traffic = "#2C5F8D", Human = "#6E9DC4", Natural = "#A9C4D8",
               Mechanical = "#C9B79C")

med_rho <- sumr$median_spearman
n_ge05  <- sumr$n_positive_0_5
write_csv(sm, file.path(DATA_DIR, "fig3_station_synchrony.csv"))

# ---- Panel a: rho distribution -----------------------------------------------
sm <- sm %>% mutate(synchrony_type = factor(synchrony_type,
                                            levels = names(SYN_COLS)))
pa <- ggplot(sm, aes(spearman_noise_heat_50m, fill = synchrony_type)) +
  geom_histogram(binwidth = 0.05, colour = "white", linewidth = 0.15,
                 boundary = 0) +
  geom_vline(xintercept = med_rho, linetype = "22", linewidth = 0.45,
             colour = "grey20") +
  scale_fill_manual(values = SYN_COLS, labels = SYN_LAB, name = NULL) +
  scale_x_continuous(breaks = seq(-0.5, 1, 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(xlim = c(-0.62, 1.0)) +
  labs(x = expression("Within-station heat-noise"~rho), y = "Stations") +
  theme_pub(BS, AT) +
  theme(legend.position = "none")

# ---- Panel b: peak-hour alignment --------------------------------------------
med_lag <- median(abs(sm$peak_lag_heat_minus_noise_hours), na.rm = TRUE)
pb <- ggplot(sm, aes(heat_peak_hour_50m, noise_peak_hour, colour = synchrony_type)) +
  annotate("segment", x = 0, y = 0, xend = 23, yend = 23,
           linetype = "22", colour = "grey55", linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.35, height = 0.35, seed = 11),
             size = 1.1, alpha = 0.8) +
  scale_colour_manual(values = SYN_COLS, labels = SYN_LAB, name = NULL) +
  scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 23.5)) +
  scale_y_continuous(breaks = seq(0, 24, 6), limits = c(0, 23.5)) +
  labs(x = "Heat-activity peak hour", y = "Noise peak hour") +
  theme_pub(BS, AT)

# ---- Panel c: high-high source composition -----------------------------------
event_map <- c("汽车行驶声" = "Traffic", "刹车声" = "Traffic", "摩托车声" = "Traffic",
               "汽车鸣笛声" = "Traffic", "交谈声" = "Human", "音乐声" = "Human",
               "叫卖声" = "Human", "店铺喇叭声" = "Human", "机械摩擦声" = "Mechanical",
               "鸟鸣声" = "Natural", "虫鸣声" = "Natural", "Unknown" = "Unknown")
comp_src <- comp %>%
  mutate(source = factor(event_map[primary_event_cn],
                         levels = c("Natural", "Mechanical", "Human", "Traffic"))) %>%
  group_by(source) %>% summarise(share = sum(share), .groups = "drop") %>%
  filter(!is.na(source))
write_csv(comp_src, file.path(DATA_DIR, "fig3_highhigh_composition.csv"))
veh_share <- comp %>% filter(primary_event_cn == "汽车行驶声") %>% pull(share)

pc <- ggplot(comp_src, aes(share, source, fill = source)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = scales::percent(share, accuracy = 1)),
            hjust = -0.22, size = 2.3, family = BASE_FAMILY, colour = "black") +
  scale_fill_manual(values = COMP_COLS, guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.93), breaks = seq(0, 0.8, 0.4),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Share of high-heat /\nhigh-noise hours", y = NULL) +
  theme_pub(BS, AT)

# ---- Compose -----------------------------------------------------------------
# Panels (a) and (b) share the synchrony-class legend -> collect once, placed
# as a single shared key on the top strip (style rule: shared legend when clean).
fig3 <- pa + pb + pc +
  plot_layout(widths = c(1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY),
        legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 6.6),
        legend.key.size = unit(0.28, "cm"), legend.margin = margin(0, 0, 0, 0))

save_fig(fig3, "fig3_synchrony.png", width_mode = "double_column", height_in = 2.7)
message("Fig 3 done. median rho=", round(med_rho,2), " median|lag|=", med_lag,
        " veh=", round(veh_share,3))
