# Exploratory figures: double-column width, square panels, shared legends.
suppressPackageStartupMessages(library(patchwork))
p11_plot_spec <- list(figure_width_mode = "double_column", fig_width_in = 178 / 25.4,
  panel_aspect = 1, font_family_primary = "Helvetica", font_family_fallback = "Arial",
  font_size_base = 9, font_size_axis_title = 10, tick_direction = "outward",
  major_tick_target = 5, minor_ticks = FALSE, show_spines = c("left", "bottom"),
  grid_mode = "none", legend_strategy = "shared top strip", panel_label_mode = "lowercase outside panel",
  output_format = "PNG", dpi = 600, transparent = FALSE)
p11_theme <- function() theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(axis.line = element_line(colour = "black", linewidth = .4),
    axis.ticks = element_line(colour = "black", linewidth = .35),
    axis.ticks.length = grid::unit(.1, "cm"), axis.title = element_text(size = 10),
    axis.text = element_text(colour = "black"), panel.grid = element_blank(),
    aspect.ratio = 1, legend.position = "top", legend.title = element_blank(),
    legend.key.width = grid::unit(9, "mm"), plot.tag = element_text(face = "bold", size = 10),
    plot.tag.position = "topleft", plot.background = element_rect(fill = "white", colour = NA))
p11_png <- function(plot, path, height_mm = 178) ggsave(path, plot,
  width = 178 / 25.4, height = height_mm / 25.4, units = "in", dpi = 600, bg = "white")
