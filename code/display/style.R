# =============================================================================
# style.R  --  Shared publication-plot style for Project 11 (metro-entrance noise)
# Follows knowledge/Academic_plot_style.md: A4 journal, single/double column,
# Helvetica 8-10 pt at print scale, left+bottom spines only, outward ticks,
# no grid, Okabe-Ito palette, PNG @ 600 dpi, opaque white background.
# Plotting library: R / ggplot2 (Project 08+ rule).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(ragg)
})

# ---- Paths (relative to the repository root; run_all.py runs every script there) ----
DATA_DIR  <- file.path("output", "figure_data")     # the per-figure CSVs each figure is drawn from
FIG_DIR   <- file.path("output", "figures")         # rendered main-text PNGs
SIFIG_DIR <- file.path("output", "figures", "si")   # rendered SI PNGs
EXPLORE_OUT   <- file.path("intermediate", "analysis")  # analysis-module result tables
EXPLORE_CLEAN <- file.path("intermediate", "clean")     # tables extracted from the deposit workbook
OSM_FILE <- file.path("data", "osm_roads_shenyang.csv") # cached OpenStreetMap road layer (Figs 1 and 5)
for (d in c(DATA_DIR, FIG_DIR, SIFIG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Width constants (mm -> in) ---------------------------------------------
WIDTH_MM  <- c(single_column = 85, double_column = 178)
mm2in     <- function(mm) as.numeric(mm) / 25.4
fig_width <- function(mode = "single_column") mm2in(WIDTH_MM[[mode]])

# ---- Typography --------------------------------------------------------------
BASE_FAMILY <- "Helvetica"   # fallback Arial handled by the system

# ---- Palettes ----------------------------------------------------------------
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7",
               "#F0E442", "#000000")

# Consistent sound-source colours used across Figs 3, 4, 5 (grayscale-safe order).
SOURCE_LEVELS <- c("Traffic", "Human", "Natural", "Mechanical", "Unknown")
SOURCE_COLS   <- c(Traffic    = "#D55E00",   # vermillion  (dominant, "hot")
                   Human      = "#0072B2",   # blue
                   Natural    = "#009E73",   # green
                   Mechanical = "#CC79A7",   # purple
                   Unknown    = "#999999")   # grey

# Diverging colour for correlation / coefficient heatmaps.
DIVERGE_LOW  <- "#0072B2"
DIVERGE_MID  <- "#FFFFFF"
DIVERGE_HIGH <- "#D55E00"

# ---- Reusable theme ----------------------------------------------------------
theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) %+replace%
    theme(
      axis.line         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.10, "cm"),                 # outward ticks
      axis.title        = element_text(size = axis_title_size, colour = "black"),
      axis.title.x      = element_text(margin = margin(t = 3)),
      axis.title.y      = element_text(margin = margin(r = 3), angle = 90),
      axis.text         = element_text(size = base_size, colour = "black"),
      legend.text       = element_text(size = base_size),
      legend.title      = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(0.34, "cm"),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      panel.grid        = element_blank(),
      plot.title        = element_text(size = axis_title_size, hjust = 0,
                                        face = "plain"),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.tag          = element_text(size = axis_title_size + 1, face = "bold",
                                       family = BASE_FAMILY),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, colour = "black")
    )
}

# Patchwork tag styling (lowercase a, b, c at top-left of full panel footprint).
tag_theme <- function(size = 11) {
  theme(plot.tag = element_text(size = size, face = "bold", family = BASE_FAMILY),
        plot.tag.position = c(0, 1))
}

# ---- Saver -------------------------------------------------------------------
# Renders PNG with ragg (good Helvetica hinting) at 600 dpi, opaque white.
save_fig <- function(plot, file, width_mode = "single_column",
                     height_in = NULL, aspect = 1, dir = FIG_DIR) {
  w <- fig_width(width_mode)
  h <- if (is.null(height_in)) w * aspect else height_in
  path <- file.path(dir, file)
  ggsave(path, plot, width = w, height = h, units = "in", dpi = 600,
         device = ragg::agg_png, bg = "white")
  message(sprintf("  wrote %s  (%.2f x %.2f in)", file, w, h))
  invisible(path)
}

# ---- Small helpers -----------------------------------------------------------
# Latitude-corrected fixed aspect ratio for a lon/lat point map (no sf needed).
map_ratio <- function(lat) 1 / cos(mean(range(lat)) * pi / 180)

# Format p-values in house style (leading zero, p < 0.001 when tiny).
fmt_p <- function(p) ifelse(p < 0.001, "italic(p) < 0.001",
                            paste0("italic(p) == ", formatC(p, digits = 3, format = "f")))

# Tidy panel-tag annotation in patchwork compositions.
PANEL_TAG <- function(p) p + theme(plot.tag = element_text(face = "bold"))

message("style.R loaded; FIG_DIR = ", FIG_DIR)
