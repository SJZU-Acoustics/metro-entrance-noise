# =============================================================================
# Figure 1 -- Study design, data and analytical framework
#   (a) Map of the 111 metro-entrance monitoring sites across Shenyang,
#       coloured by metro line; entrances without valid 24 h acoustic data
#       shown as open symbols.
#   (b) Data architecture: the five measured/derived domains feeding the
#       station-hour record and the four analytical lines of the paper.
#   (c) Analytical framework (Kang reframe 2026-07-01; mirrored layout per
#       Prof Zhang's sketch 2026-07-06): noise at the city's transit activity
#       nodes (RQ1, outcome box on the left), driven in time by human activity
#       (RQ2, from above), explained by sound-source composition vs the
#       spatial-physical context of the node setting (RQ3), and the context-
#       composition relationship examined in RQ4. Context branches
#       follow the soundscape-context framework (Zhang et al. 2025, JASA).
#       Box fills reuse panel (b)'s domain colours (noise = sound levels,
#       activity = population heat, context = built environment / greenness);
#       the sound-source composition box is left white (the construct set
#       against context in RQ3/RQ4), and the out-of-scope socio-cultural branch is a
#       grey box.
# =============================================================================
source(file.path("code", "display", "style.R"))

# ---- Data --------------------------------------------------------------------
sm <- read_csv(file.path(EXPLORE_CLEAN, "station_master.csv"), show_col_types = FALSE)
nm <- read_csv(file.path(EXPLORE_CLEAN, "noise_metrics_station.csv"), show_col_types = FALSE)

line_label <- function(x) {
  num <- str_match(x, "地铁(\\d+)号线")[, 2]
  factor(paste0("Line ", num), levels = paste0("Line ", c(1, 2, 4, 9, 10)))
}

map_df <- sm %>%
  transmute(station_key, lng = lng_round, lat = lat_round,
            line = line_label(line_name)) %>%
  left_join(nm %>% select(station_key, laeq_24h_leq_db_a), by = "station_key") %>%
  mutate(acoustic = ifelse(is.na(laeq_24h_leq_db_a), "No valid acoustic data",
                           "Valid 24 h record"))
write_csv(map_df, file.path(DATA_DIR, "fig1_station_map.csv"))

n_valid   <- sum(map_df$acoustic == "Valid 24 h record")
n_missing <- sum(map_df$acoustic != "Valid 24 h record")

LINE_COLS <- c("Line 1" = "#0072B2", "Line 2" = "#D55E00", "Line 4" = "#009E73",
               "Line 9" = "#CC79A7", "Line 10" = "#E69F00")

# ---- Panel a: entrance map ---------------------------------------------------
rr <- map_ratio(map_df$lat)
missing_df <- map_df %>% filter(acoustic == "No valid acoustic data")
xr <- range(map_df$lng) + c(-0.015, 0.015); yr <- range(map_df$lat) + c(-0.012, 0.012)
# Major-road context from the cached OpenStreetMap layer (code/fetch_osm.py);
# drawn only when the cache exists.
road_file <- OSM_FILE
roads <- if (file.exists(road_file)) {
  read_csv(road_file, show_col_types = FALSE) %>%
    filter(lng > xr[1] - 0.02, lng < xr[2] + 0.02, lat > yr[1] - 0.02, lat < yr[2] + 0.02)
} else NULL
road_layer <- if (is.null(roads)) NULL else
  geom_path(data = roads, aes(lng, lat, group = line_id), colour = "grey82",
            linewidth = 0.25, inherit.aes = FALSE)
pa <- ggplot(map_df, aes(lng, lat)) +
  road_layer +
  geom_point(aes(colour = line), size = 1.7, stroke = 0.6) +
  geom_point(data = missing_df, aes(lng, lat),
             shape = 1, size = 2.6, stroke = 0.7, colour = "grey25") +
  ggrepel::geom_text_repel(data = missing_df,
                           aes(lng, lat, label = station_key),
                           size = 2.1, family = BASE_FAMILY, colour = "grey20",
                           segment.size = 0.3, segment.colour = "grey55",
                           min.segment.length = 0, force = 1.2,
                           direction = "both", max.overlaps = Inf, seed = 11) +
  scale_colour_manual(values = LINE_COLS, name = NULL) +
  scale_x_continuous(breaks = seq(123.25, 123.50, 0.10)) +
  scale_y_continuous(breaks = seq(41.65, 41.95, 0.10)) +
  coord_fixed(ratio = rr, xlim = xr, ylim = yr, expand = FALSE) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") +
  theme_pub(base_size = 7.5, axis_title_size = 8.5) +
  theme(legend.position = "right",
        legend.key.size = unit(0.30, "cm"),
        legend.spacing.y = unit(0.02, "cm"),
        legend.text = element_text(size = 6.6),
        legend.margin = margin(0, 0, 0, -4))

# ---- Panel b: data-architecture schematic ------------------------------------
# Native-R reproduction (2026-07-02, for reproducibility) of Prof Zhang's
# supplied raster redesign (assets/fig1b_data_architecture.png, kept as the
# visual reference). Geometry, fills, text colours and line weights were
# measured from that PNG (1030 x 1246 px) and are drawn here in its pixel
# coordinates (y increases downward via scale_y_reverse; coord_fixed keeps the
# native aspect, so patchwork letterboxes it exactly as the raster was).
BB_GREY <- "#595959"   # borders, arrows, body text (glyph colour 89,89,89)
bb_box <- function(x0, x1, y0, y1, fill)
  annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1,
           fill = fill, colour = BB_GREY, linewidth = 0.33)
bb_txt <- function(x, y, label, size = 2.45, face = "plain", hjust = 0.5,
                   col = BB_GREY, angle = 0, parse = FALSE)
  annotate("text", x = x, y = y, label = label, size = size, family = BASE_FAMILY,
           fontface = face, hjust = hjust, colour = col, angle = angle,
           lineheight = 0.95, parse = parse)
bb_arr <- function(x1, y1, x2, y2, lw = 0.51)
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2,
           arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
           linewidth = lw, colour = BB_GREY, linejoin = "mitre")

# The five domain boxes: (top, bottom) in pixel rows; fills sampled from the PNG.
bb_dom <- tibble(
  y0   = c(38, 227, 421, 610, 804),
  y1   = c(214, 408, 597, 791, 923),
  fill = c("#DAE3F3", "#EDF6FC", "#EDE4EA", "#D4E6E7", "#D4E6E7"),
  head = c("Sound levels", "Sound events", "Population heat",
           "Built environment", "Greenness"))

pb <- ggplot() +
  scale_x_continuous(limits = c(0, 1030), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1246, 0), expand = c(0, 0)) +
  # -- dashed rounded container --------------------------------------------------
  annotation_custom(grid::roundrectGrob(
    r = unit(2.6, "mm"),
    gp = grid::gpar(col = BB_GREY, fill = NA, lwd = 1.0, lty = "dashed")),
    xmin = 322, xmax = 1017, ymin = -942, ymax = -18) +
  # -- left source box (rotated text) ---------------------------------------------
  bb_box(28, 232, -942, -18, fill = "#F2F2F2") +
  bb_txt(94, -480, "109 metro entrances", size = 2.55, face = "bold",
         col = "black", angle = 90) +
  bb_txt(175, -500, "24 h field monitoring; 2,616 station-hours",
         size = 2.4, angle = 90) +
  # -- five domain boxes -----------------------------------------------------------
  { dl <- list(); for (i in seq_len(nrow(bb_dom))) {
      dl <- c(dl, list(
        bb_box(383, 963, -bb_dom$y1[i], -bb_dom$y0[i], fill = bb_dom$fill[i]),
        bb_txt(417, -(bb_dom$y0[i] + 33), bb_dom$head[i], size = 2.55,
               face = "bold", col = "black", hjust = 0)))
    }; dl } +
  # body lines (centred on x = 673; y at measured band midpoints)
  bb_txt(673, -134, "italic(L)['Aeq']*' (24h, day, night);'", parse = TRUE) +
  bb_txt(673, -188, "diurnal amplitude") +
  bb_txt(673, -317, "Traffic / human / nature /") +
  bb_txt(673, -373, "mechanical share") +
  bb_txt(673, -515, "50m hourly activity") +
  bb_txt(673, -571, "proxy (24h)") +
  bb_txt(673, -704, "Land use, road density,") +
  bb_txt(673, -760, "POI (100/300/500m)") +
  bb_txt(673, -898, "NDVI (20-500m)") +
  # -- fan arrows: left box -> each domain box --------------------------------------
  { al <- list(); for (i in seq_len(nrow(bb_dom)))
      al <- c(al, list(bb_arr(238, -480, 376, -(bb_dom$y0[i] + bb_dom$y1[i])/2)));
    al } +
  # -- container -> station-level analysis ------------------------------------------
  bb_arr(670, -946, 670, -1003, lw = 0.62) +
  bb_box(36, 1017, -1231, -1011, fill = "#F2F2F2") +
  bb_txt(526, -1060, "Station-level analysis", size = 2.55, face = "bold",
         col = "black") +
  bb_txt(526, -1109, "Temporal burden · activity synchrony ·") +
  bb_txt(526, -1175, "source mechanism · spatial context") +
  coord_fixed(ratio = 1, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 10))

# (raster reference version, superseded 2026-07-02)
if (FALSE) {
img_b <- png::readPNG(file.path(P11_ROOT, "assets", "fig1b_data_architecture.png"))
pb <- ggplot() +
  annotation_custom(grid::rasterGrob(img_b, interpolate = TRUE)) +
  theme_void() +
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 10))
}

# ---- (deprecated) hand-laid schematic on a 0-100 canvas ----------------------
if (FALSE) {
box <- function(x, y, w, h, fill, col = "grey30", r = NA)
  annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
           fill = fill, colour = col, linewidth = 0.4)
txt <- function(x, y, label, size = 2.5, face = "plain", hjust = 0.5, col = "black")
  annotate("text", x = x, y = y, label = label, size = size, family = BASE_FAMILY,
           fontface = face, hjust = hjust, colour = col, lineheight = 0.95)
arr <- function(x1, y1, x2, y2)
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2,
           arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
           linewidth = 0.4, colour = "grey35")

dom_y <- seq(88, 22, length.out = 5)
dom <- tibble(
  y    = dom_y,
  fill = c("#FBE6D4", "#E7EBF3", "#D9EFE6", "#F3E2EE", "#EDEFF0"),
  head = c("Acoustic levels", "Sound events", "Population heat",
           "Built environment", "Greenness"),
  body = c("LAeq 24 h, day, night;\ndiurnal amplitude",
           "traffic / human / natural /\nmechanical share",
           "50 m hourly activity\nproxy (24 h)",
           "land use, road density,\nPOI (100/300/500 m)",
           "NDVI (20-500 m)"))

# Source node: centred vertically on the domain stack.
src_x <- 14; src_y <- 55
# Domain boxes stacked on the right, spanning the vertical extent.
dom_x <- 60; dom_w <- 42; dom_h <- 11.5
# Analysis node: full-width band at the bottom.
an_y <- 6.5

pb <- ggplot() +
  xlim(0, 100) + ylim(0, 100) +
  # -- left node: monitoring campaign ----------------------------------------
  box(src_x, src_y, 22, 34, fill = "#F2F2F2", col = "grey25") +
  txt(src_x, src_y + 11.5, "111 metro", size = 2.9, face = "bold") +
  txt(src_x, src_y + 7.5,  "entrances", size = 2.9, face = "bold") +
  txt(src_x, src_y + 1,    "24 h field\nmonitoring", size = 2.5) +
  txt(src_x, src_y - 7.5,  "2,664\nstation-hours", size = 2.3, col = "grey30") +
  # -- five domain boxes -----------------------------------------------------
  { dl <- list(); for (i in seq_len(nrow(dom))) {
      dl <- c(dl,
        box(dom_x, dom$y[i], dom_w, dom_h, fill = dom$fill[i]),
        txt(dom_x - dom_w/2 + 2.6, dom$y[i] + 2.6, dom$head[i],
            size = 2.55, face = "bold", hjust = 0),
        txt(dom_x - dom_w/2 + 2.6, dom$y[i] - 2.7, dom$body[i],
            size = 2.1, hjust = 0, col = "grey25"))
    }; dl } +
  # -- bottom node: four analytical lines ------------------------------------
  box(dom_x, an_y, 78, 12, fill = "#F2F2F2", col = "grey25") +
  txt(dom_x, an_y + 3.0, "Station-level analysis", size = 2.7, face = "bold") +
  txt(dom_x, an_y - 2.4,
      "temporal burden  ·  activity synchrony  ·  source mechanism  ·  spatial context",
      size = 2.1, col = "grey25") +
  # -- arrows: monitoring -> each domain -------------------------------------
  { al <- list(); for (i in seq_len(nrow(dom)))
      al <- c(al, arr(src_x + 11, src_y, dom_x - dom_w/2, dom$y[i])); al } +
  # -- arrows: domains converge -> analysis ----------------------------------
  { bl <- list(); for (i in seq_len(nrow(dom)))
      bl <- c(bl, arr(dom_x, dom$y[i] - dom_h/2, dom_x, an_y + 6)); bl } +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 8) +
  theme(plot.margin = margin(2, 2, 2, 2))
}

# ---- Panel c: analytical framework -------------------------------------------
# Mirrored layout (Prof Zhang's PowerPoint sketch, 2026-07-06) on the same
# 0-100 x -1..41 canvas, styled to match panel (b): the outcome box sits on
# the left, driven from above by human activity (RQ2) and from the right by
# the two compared explanations (RQ3) -- sound-source composition (upper path)
# and the context of the node setting -- with the context-composition
# relationship examined in RQ4. Fills reuse panel (b)'s domain colours so
# each construct matches the data domain that measures it, except the central
# sound-source composition box, left white to set it apart as the construct
# set against context in RQ3/RQ4; the out-of-scope socio-cultural box is grey like the
# neutral source/analysis boxes.
fbox <- function(x, y, w, h, fill, col = "grey30", lty = "solid")
  annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
           fill = fill, colour = col, linewidth = 0.4, linetype = lty)
ftxt <- function(x, y, label, size = 2.2, face = "plain", hjust = 0.5, col = "grey10")
  annotate("text", x = x, y = y, label = label, size = size, family = BASE_FAMILY,
           fontface = face, hjust = hjust, colour = col, lineheight = 1.0)
farr <- function(x1, y1, x2, y2)
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2,
           arrow = arrow(length = unit(0.19, "cm"), type = "closed"),
           linewidth = 0.55, colour = "grey30", linejoin = "mitre")

pc <- ggplot() +
  xlim(0, 100) + ylim(-1, 44) +
  # -- human activity (top left) --------------------------------------------------
  fbox(16, 40, 28, 7, fill = "#EDE4EA") +
  ftxt(16, 41.45, "Human activity", size = 2.35, face = "bold") +
  ftxt(16, 38.75, "hourly population presence", size = 1.95, col = "grey35") +
  # -- noise at the node: the outcome (left) ---------------------------------------
  fbox(16, 19, 28, 26, fill = "#DAE3F3") +
  ftxt(16, 28.4, "Noise at the node  (RQ1)", size = 2.55, face = "bold") +
  ftxt(4.4, 24.9, "high · time-structured · uneven", size = 2.1, hjust = 0) +
  annotate("text", x = 4.4, y = 20.3, hjust = 0, size = 2.1, family = BASE_FAMILY,
           colour = "grey10",
           label = "'level  '*italic(L)['Aeq']*' (24h, day, night)'",
           parse = TRUE) +
  ftxt(4.4, 16.2, "diurnal cycle · amplitude", size = 2.1, hjust = 0) +
  ftxt(4.4, 12.1, "spatial clustering", size = 2.1, hjust = 0) +
  # -- sound-source composition (middle): white -- the pivotal mediator ------------
  fbox(51, 25, 28, 10, fill = "white") +
  ftxt(51, 27.3, "Sound-source composition", size = 2.55, face = "bold") +
  ftxt(51, 23.0, "traffic · human · natural · mechanical", size = 2.05, col = "grey25") +
  # -- context of the node setting (right) -----------------------------------------
  fbox(85, 19, 28, 26, fill = "#D4E6E7") +
  ftxt(85, 28.4, "Context of the node setting", size = 2.55, face = "bold") +
  ftxt(85, 24.9, "spatial–physical attributes", size = 1.95, col = "grey35") +
  ftxt(73.4, 20.3, "Macro – land use · road network", size = 2.1, hjust = 0) +
  ftxt(73.4, 16.2, "Meso – greenness · facilities", size = 2.1, hjust = 0) +
  ftxt(73.4, 12.1, "Node – metro line · station type", size = 2.1, hjust = 0) +
  # -- socio-cultural branch: out of scope, grey (top right) ------------------------
  fbox(85, 40, 28, 7, fill = "#F2F2F2", col = "grey45") +
  ftxt(85, 41.45, "Socio-cultural attributes", size = 2.35, face = "bold",
       col = "grey30") +
  ftxt(85, 38.75, "shape perception, not measured here", size = 1.95, col = "grey45") +
  # -- arrows (one grey family, as in panel b) --------------------------------------
  farr(16, 36.5, 16, 32.4) +                      # activity -> noise (RQ2)
  ftxt(14.4, 34.2, "RQ2", size = 2.3, face = "bold", hjust = 1, col = "grey15") +
  farr(37, 25, 30.8, 25) +                        # composition -> noise (RQ3)
  ftxt(33.9, 27.2, "RQ3", size = 2.3, face = "bold", col = "grey15") +
  farr(71, 25, 65.8, 25) +                        # context -> composition (RQ4)
  ftxt(68.4, 27.2, "RQ4", size = 2.3, face = "bold", col = "grey15") +
  farr(71, 9.5, 30.8, 9.5) +                      # context -> noise (RQ3)
  ftxt(50.9, 11.5, "RQ3", size = 2.3, face = "bold", col = "grey15") +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 8) +
  theme(plot.margin = margin(4, 2, 2, 2))

# ---- Compose -----------------------------------------------------------------
fig1 <- (pa + pb + plot_layout(widths = c(1.5, 1))) / pc +
  plot_layout(heights = c(2.1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY))

save_fig(fig1, "fig1_study_map.png", width_mode = "double_column", height_in = 5.1)
message("Fig 1 done: ", n_valid, " valid / ", n_missing, " missing acoustic.")
