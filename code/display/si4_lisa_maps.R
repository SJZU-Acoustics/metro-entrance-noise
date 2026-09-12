# SI Fig S3 -- Local spatial clustering (LISA) of entrance noise.
# Two labelled panels (daytime, night-time); shared cluster legend.
source(file.path("code", "display", "style.R"))

sm <- read_csv(file.path(EXPLORE_CLEAN, "station_master.csv"), show_col_types = FALSE) %>%
  select(station_key, lng = lng_round, lat = lat_round)
read_lisa <- function(file, lab) {
  read_csv(file.path(EXPLORE_OUT, file), show_col_types = FALSE) %>%
    transmute(station_key, cluster = lisa_cluster, period = lab) %>%
    left_join(sm, by = "station_key")
}
li <- bind_rows(
  read_lisa("12_spatial__lisa__day__lday_leq_db_a.csv", "Daytime (Lday)"),
  read_lisa("12_spatial__lisa__night__lnight_leq_db_a.csv", "Night-time (Lnight)")
) %>% mutate(cluster = factor(cluster, levels = c("HH","LL","HL","LH","NS")))

BS <- 8; AT <- 9
CL_COLS <- c(HH = "#D55E00", LL = "#0072B2", HL = "#E69F00", LH = "#56B4E9", NS = "grey80")
CL_LAB  <- c(HH = "High-High", LL = "Low-Low", HL = "High-Low", LH = "Low-High",
             NS = "Not significant")
rr <- map_ratio(li$lat)

map_panel <- function(per, show_y) {
  p <- ggplot(filter(li, period == per), aes(lng, lat)) +
    geom_point(aes(colour = cluster, size = cluster == "NS")) +
    scale_colour_manual(values = CL_COLS, labels = CL_LAB, name = NULL) +
    scale_size_manual(values = c(`TRUE` = 0.8, `FALSE` = 1.7), guide = "none") +
    scale_x_continuous(breaks = seq(123.25, 123.50, 0.1)) +
    scale_y_continuous(breaks = seq(41.65, 41.95, 0.1)) +
    coord_fixed(ratio = rr) +
    guides(colour = guide_legend(override.aes = list(size = 1.8))) +
    labs(x = "Longitude (°E)", y = if (show_y) "Latitude (°N)" else NULL) +
    theme_pub(BS, AT)
  if (!show_y) p <- p + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  p
}
pa <- map_panel("Daytime (Lday)", TRUE)
pb <- map_panel("Night-time (Lnight)", FALSE)

fig <- pa + pb +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY),
        legend.position = "right", legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))
save_fig(fig, "si_lisa_maps.png", width_mode = "double_column", height_in = 3.1, dir = SIFIG_DIR)
message("SI S3 done.")
