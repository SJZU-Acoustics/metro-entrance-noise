# SI Fig S2 -- Multiscale built-environment partial-Spearman associations.
# One labelled panel per outcome (no facet-strip titles); shared colour legend.
source(file.path("code", "display", "style.R"))

d <- read_csv(file.path(EXPLORE_OUT, "06_built_env__correlations.csv"), show_col_types = FALSE)
BS <- 8; AT <- 9

cat_lab <- c(road_arterial_density = "Arterial road density",
             poi_total_density = "POI density",
             commercial_land_pct = "Commercial land",
             landuse_mix_shannon = "Land-use mix",
             residential_land_pct = "Residential land",
             road_pedestrian_density = "Pedestrian road density",
             road_expressway_density = "Expressway density",
             green_land_pct = "Green-space land",
             water_body_pct = "Water body")
out_lab <- c(laeq_24h_leq_db_a = "LAeq 24 h", lday_leq_db_a = "Lday",
             lnight_leq_db_a = "Lnight")

hm <- d %>%
  filter(method == "partial_spearman", outcome %in% names(out_lab)) %>%
  mutate(category = factor(cat_lab[category], levels = rev(unname(cat_lab))),
         outcome = out_lab[outcome],
         scale = factor(paste0(scale_m, " m"), levels = c("100 m","300 m","500 m")),
         sig = q_fdr < 0.05) %>%
  filter(!is.na(category))

heat_panel <- function(oc) {
  ggplot(filter(hm, outcome == oc), aes(scale, category, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_point(data = ~filter(.x, sig), shape = 21, size = 0.9, stroke = 0.5,
               colour = "black", fill = "black") +
    scale_fill_gradient2(low = DIVERGE_LOW, mid = DIVERGE_MID, high = DIVERGE_HIGH,
                         midpoint = 0, limits = c(-0.4, 0.4), name = "Partial ρ") +
    labs(x = "Buffer radius", y = NULL) +
    theme_pub(BS, AT) +
    theme(axis.line = element_blank(), axis.ticks = element_line(linewidth = 0.3),
          axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7))
}
pa <- heat_panel("LAeq 24 h")
pb <- heat_panel("Lday")
pc <- heat_panel("Lnight")

fig <- pa + pb + pc +
  plot_layout(guides = "collect", axes = "collect", widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY),
        legend.position = "right", legend.title = element_text(size = 7.5, angle = 90),
        legend.title.position = "left")
save_fig(fig, "si_builtenv_multiscale.png", width_mode = "double_column", height_in = 2.9, dir = SIFIG_DIR)
message("SI S2 done.")
