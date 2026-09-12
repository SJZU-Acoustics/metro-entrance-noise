# =============================================================================
# Figure 5 -- Where node noise and natural sound sit, the spatial models,
#             and the greenness lever
#   (a) Point map of LAeq_24h at the 109 entrances (two invalid entrances open).
#   (b) Point map of natural-source share on the same extent.
#   (c) Source-term coefficients (95% CI) for LAeq_24h across OLS,
#       trend-surface, spatial-lag and spatial-error models.
#   (d) Residual Moran's I for LAeq_24h as a source share and a spatial trend
#       are added.
#   (e) Natural-source share vs NDVI (500 m), partial correlation given LAeq_24h.
# =============================================================================
source(file.path("code", "display", "style.R"))

sm <- read_csv(file.path(EXPLORE_CLEAN, "station_master.csv"), show_col_types = FALSE)
nm <- read_csv(file.path(EXPLORE_CLEAN, "noise_metrics_station.csv"), show_col_types = FALSE)
coef <- read_csv(file.path(EXPLORE_OUT, "28_source_spatial__coefficients.csv"),
                 show_col_types = FALSE)
moran <- read_csv(file.path(EXPLORE_OUT, "28_source_spatial__residual_moran.csv"),
                  show_col_types = FALSE)
div  <- read_csv(file.path(EXPLORE_OUT, "24_acoustic_diversity__station_metrics.csv"),
                 show_col_types = FALSE)
boot <- read_csv(file.path(EXPLORE_OUT, "29_beyond_loudness_partial__boot_ci.csv"),
                 show_col_types = FALSE)

BS <- 8; AT <- 9
TERM_COLS <- c("Traffic share" = "#D55E00", "Natural-source share" = "#009E73")

# ---- Map data ----------------------------------------------------------------
map_df <- sm %>%
  transmute(station_key, lng = lng_round, lat = lat_round, natural = sound_natural_share) %>%
  left_join(nm %>% select(station_key, laeq = laeq_24h_leq_db_a), by = "station_key")
write_csv(map_df, file.path(DATA_DIR, "fig5_station_maps.csv"))
valid   <- map_df %>% filter(!is.na(laeq))
missing <- map_df %>% filter(is.na(laeq))
rr <- map_ratio(map_df$lat)

# Optional road context: a cached OpenStreetMap arterial-road layer
# (data/osm_roads_shenyang.csv: line_id, lng, lat, drawn as grey paths).
road_file <- OSM_FILE
xr <- range(map_df$lng) + c(-0.015, 0.015); yr <- range(map_df$lat) + c(-0.012, 0.012)
roads <- if (file.exists(road_file)) {
  read_csv(road_file, show_col_types = FALSE) %>%
    filter(lng > xr[1] - 0.02, lng < xr[2] + 0.02, lat > yr[1] - 0.02, lat < yr[2] + 0.02)
} else NULL
road_layer <- function() {
  if (is.null(roads)) return(NULL)
  geom_path(data = roads, aes(lng, lat, group = line_id), colour = "grey82",
            linewidth = 0.25, inherit.aes = FALSE)
}

base_map <- function(d, aes_col, legend_name, low, high, labels = waiver()) {
  ggplot(d, aes(lng, lat)) +
    road_layer() +
    geom_point(aes(colour = {{ aes_col }}), size = 1.9) +
    geom_point(data = missing, aes(lng, lat), shape = 1, size = 2.4, stroke = 0.6,
               colour = "grey35") +
    scale_colour_gradient(low = low, high = high, name = legend_name, labels = labels) +
    scale_x_continuous(breaks = seq(123.25, 123.50, 0.10)) +
    scale_y_continuous(breaks = seq(41.65, 41.95, 0.10)) +
    coord_fixed(ratio = rr, xlim = xr, ylim = yr, expand = FALSE) +
    labs(x = "Longitude (°E)", y = "Latitude (°N)") +
    theme_pub(BS, AT) +
    theme(legend.position = "right",
          legend.title = element_text(size = 7.2),
          legend.text = element_text(size = 6.8),
          legend.key.height = unit(0.42, "cm"), legend.key.width = unit(0.28, "cm"),
          legend.margin = margin(0, 0, 0, 0))
}
pa <- base_map(valid, laeq, expression(italic(L)[plain("Aeq,24h")]~"(dB(A))"),
               low = "#56B4E9", high = "#D55E00")
pb <- base_map(valid, natural, "Natural-source share", low = "grey88", high = "#009E73",
               labels = scales::percent_format(accuracy = 1))

# ---- Panel c: coefficient forest, LAeq_24h ----------------------------------
MODEL_LAB <- c(OLS = "OLS", OLS_trend = "Trend surface",
               SAR = "Spatial lag", SEM = "Spatial error")
fr <- coef %>%
  filter(outcome == "LAeq_24h",
         term %in% c("sound_traffic_share", "sound_natural_share")) %>%
  transmute(outcome, model,
            term = recode(term, sound_traffic_share = "Traffic share",
                          sound_natural_share = "Natural-source share"),
            coef, lo = coef - 1.96 * se, hi = coef + 1.96 * se) %>%
  mutate(model = factor(MODEL_LAB[model], levels = rev(MODEL_LAB)),
         term  = factor(term, levels = names(TERM_COLS)))
write_csv(fr, file.path(DATA_DIR, "fig5_coef_forest.csv"))
pc <- ggplot(fr, aes(coef, model, colour = term)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey55", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(size = 1.3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = TERM_COLS, name = NULL) +
  scale_x_continuous(breaks = seq(-20, 10, 10), limits = c(-21, 14)) +
  labs(x = "Coefficient (dB per unit share)", y = NULL) +
  theme_pub(BS, AT) +
  theme(legend.position = "none")   # colour legend shared with panel d

# ---- Panel d: residual Moran's I, LAeq_24h ----------------------------------
STEP_LV <- c("Intercept\nonly", "+ source", "+ source\n+ trend")
mo <- moran %>%
  filter(outcome == "LAeq_24h") %>%
  mutate(src  = case_when(spec == "intercept_only" ~ "both",
                          grepl("traffic", spec) ~ "Traffic share",
                          TRUE ~ "Natural-source share"),
         step = case_when(spec == "intercept_only" ~ STEP_LV[1],
                          grepl("trend", spec) ~ STEP_LV[3],
                          TRUE ~ STEP_LV[2]))
mo <- bind_rows(mo %>% filter(src != "both"),
                mo %>% filter(src == "both") %>% mutate(src = "Traffic share"),
                mo %>% filter(src == "both") %>% mutate(src = "Natural-source share")) %>%
  mutate(step = factor(step, levels = STEP_LV), src = factor(src, levels = names(TERM_COLS)))
write_csv(mo, file.path(DATA_DIR, "fig5_resid_moran.csv"))
pd <- ggplot(mo, aes(step, resid_moran_I, colour = src, group = src)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = TERM_COLS, name = NULL) +
  scale_x_discrete(expand = expansion(add = 0.45)) +
  scale_y_continuous(breaks = seq(0, 0.25, 0.1), limits = c(-0.04, 0.27)) +
  labs(x = NULL, y = "Residual Moran's I") +
  theme_pub(BS, AT) +
  theme(axis.text.x = element_text(size = 6.4, lineheight = 0.9),
        legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1), legend.text = element_text(size = 6.6),
        legend.key.size = unit(0.28, "cm"))

# ---- Panel e: greenness lever ------------------------------------------------
pcb <- boot %>% filter(association == "sound_natural_share ~ ndvi_500 | laeq_24h_leq_db_a")
pe <- ggplot(div, aes(ndvi_500, sound_natural_share)) +
  geom_point(colour = "#009E73", size = 1.2, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey20",
              fill = "grey75", linewidth = 0.5, alpha = 0.35) +
  annotate("text", x = 0.12, y = 0.97, hjust = 0, vjust = 1, parse = TRUE,
           label = sprintf("atop('partial '*rho*' | '*italic(L)[plain('Aeq,24h')], '= %.2f [%.2f, %.2f]')",
                           pcb$rho, pcb$ci_lo, pcb$ci_hi),
           size = 2.1, family = BASE_FAMILY, colour = "grey20", lineheight = 0.95) +
  scale_x_continuous(breaks = seq(0.1, 0.5, 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "NDVI (500 m)", y = "Natural-source share") +
  theme_pub(BS, AT)

# ---- Compose: two maps above three analytical panels -------------------------
top    <- pa + pb + plot_layout(ncol = 2)
bottom <- pc + pd + pe + plot_layout(ncol = 3, widths = c(1.05, 1, 1))
fig5 <- (top / bottom) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY))

save_fig(fig5, "fig5_spatial_greenness.png", width_mode = "double_column", height_in = 5.6)
message("Fig 5 done (two maps + three panels).")
