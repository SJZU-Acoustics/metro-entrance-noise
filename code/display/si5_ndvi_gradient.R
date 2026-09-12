# SI Fig S6 -- Greenness (NDVI) gradient of the entrance soundscape by buffer scale.
source(file.path("code", "display", "style.R"))

co  <- read_csv(file.path(EXPLORE_OUT, "25_ndvi_gradient__correlations.csv"),
                show_col_types = FALSE)
div <- read_csv(file.path(EXPLORE_OUT, "24_acoustic_diversity__station_metrics.csv"),
                show_col_types = FALSE)
BS <- 8; AT <- 9

# Panel a: rho by buffer radius for three targets
keep <- c(laeq_24h_leq_db_a = "LAeq 24 h", lnight_leq_db_a = "Lnight",
          sound_natural_share = "Natural-source share")
ga <- co %>% filter(outcome %in% names(keep)) %>%
  mutate(target = factor(keep[outcome], levels = unname(keep)),
         sig = q_fdr < 0.05)
GCOL <- c("LAeq 24 h" = "#0072B2", "Lnight" = "#CC79A7",
          "Natural-source share" = "#009E73")
pa <- ggplot(ga, aes(ndvi_scale_m, rho, colour = target, group = target)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(shape = sig), size = 1.5) +
  scale_colour_manual(values = GCOL, name = NULL,
                      labels = c(expression(italic(L)[plain("Aeq")]~"24 h"),
                                 expression(italic(L)[plain("night")]),
                                 expression("Natural-source share"))) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_x_continuous(breaks = c(20, 100, 200, 300, 500)) +
  labs(x = "NDVI buffer radius (m)", y = expression("Spearman"~rho)) +
  theme_pub(BS, AT) +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 6.8), legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0, 0, 0, 0))

# Panel b: NDVI 500 m vs LAeq_24h (the quieter-greener side)
rho_lq <- co %>% filter(ndvi_scale_m == 500, outcome == "laeq_24h_leq_db_a") %>% pull(rho)
pb <- ggplot(div, aes(ndvi_500, laeq_24h_leq_db_a)) +
  geom_point(colour = "#0072B2", size = 1.2, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey20",
              fill = "grey75", linewidth = 0.5, alpha = 0.35) +
  annotate("text", x = 0.12, y = 52, hjust = 0,
           label = sprintf("ρ = %.2f (500 m)", rho_lq),
           size = 2.3, family = BASE_FAMILY, colour = "grey20") +
  scale_x_continuous(breaks = seq(0.1, 0.5, 0.1)) +
  scale_y_continuous(breaks = seq(50, 75, 5)) +
  labs(x = "NDVI (500 m)", y = expression(italic(L)[plain("Aeq")]~"24 h (dB(A))")) +
  theme_pub(BS, AT)

fig <- pa + pb + plot_layout(widths = c(1.05, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY))
save_fig(fig, "si_ndvi_gradient.png", width_mode = "double_column", height_in = 2.9, dir = SIFIG_DIR)
message("SI S6 done.")
