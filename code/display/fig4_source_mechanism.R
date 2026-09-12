# =============================================================================
# Figure 4 -- Sound-source composition as the dominant explanation
#   (a) Distribution of station source shares (traffic/human/natural/mechanical).
#   (b) Traffic share vs LAeq_24h.
#   (c) Natural-source share vs LAeq_24h.
#   (d) Repeated 5-fold cross-validated RMSE by predictor domain and outcome.
# =============================================================================
source(file.path("code", "display", "style.R"))

div <- read_csv(file.path(EXPLORE_OUT, "24_acoustic_diversity__station_metrics.csv"),
                show_col_types = FALSE)
shares <- read_csv(file.path(EXPLORE_CLEAN, "station_master.csv"), show_col_types = FALSE)
corr <- read_csv(file.path(EXPLORE_OUT, "08_sound_source__corr_spearman.csv"),
                 show_col_types = FALSE)
cv  <- read_csv(file.path(EXPLORE_OUT, "21_domain_cv__model_performance.csv"),
                show_col_types = FALSE)
domc <- read_csv(file.path(EXPLORE_OUT, "08_sound_source__dominant_event_counts.csv"),
                 show_col_types = FALSE)

BS <- 8.5; AT <- 9.5
get_rho <- function(out, pred) corr %>%
  filter(outcome == out, predictor == pred) %>% slice(1) %>% pull(spearman_rho)

# ---- Panel a: source-share distribution --------------------------------------
shr <- shares %>%
  filter(station_key %in% div$station_key) %>%
  transmute(station_key,
            Traffic = sound_traffic_share, Human = sound_human_share,
            Natural = sound_natural_share, Mechanical = sound_mechanical_share) %>%
  pivot_longer(-station_key, names_to = "source", values_to = "share") %>%
  mutate(source = factor(source, levels = c("Mechanical", "Human", "Natural", "Traffic")))
write_csv(shr, file.path(DATA_DIR, "fig4_source_shares.csv"))
n_veh <- domc %>% filter(dominant_event_en == "Vehicle pass-by") %>% pull(station_count)

pa <- ggplot(shr, aes(share, source, fill = source)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.35, width = 0.62,
               outlier.alpha = 0.5) +
  scale_fill_manual(values = SOURCE_COLS, guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  labs(x = "Share of hours by dominant source", y = NULL) +
  theme_pub(BS, AT)

# ---- Panels b, c: share vs LAeq_24h ------------------------------------------
scatter_panel <- function(xvar, xlab, src_col, rho, ann_extra = NULL) {
  rlab <- sprintf("ρ = %.2f, p < 0.001", rho)
  p <- ggplot(div, aes(.data[[xvar]], laeq_24h_leq_db_a)) +
    geom_point(colour = src_col, size = 1.3, alpha = 0.75) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "grey20",
                fill = "grey75", linewidth = 0.5, alpha = 0.35) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(breaks = seq(50, 75, 5)) +
    labs(x = xlab, y = expression(italic(L)[plain("Aeq")]~"24 h (dB(A))")) +
    theme_pub(BS, AT)
  list(p = p, rlab = rlab)
}
sb <- scatter_panel("sound_traffic_share", "Traffic share", SOURCE_COLS["Traffic"],
                    get_rho("laeq_24h_leq_db_a", "sound_traffic_share"))
pb <- sb$p + annotate("text", x = 0.02, y = 73.5, hjust = 0, label = sb$rlab,
                      size = 2.3, family = BASE_FAMILY, colour = "grey15")

rho_nat_night <- get_rho("lnight_leq_db_a", "sound_natural_share")
sc <- scatter_panel("sound_natural_share", "Natural-source share", SOURCE_COLS["Natural"],
                    get_rho("laeq_24h_leq_db_a", "sound_natural_share"))
pc <- sc$p +
  annotate("text", x = 0.98, y = 73.5, hjust = 1, label = sc$rlab,
           size = 2.3, family = BASE_FAMILY, colour = "grey15")

# ---- Panel d: domain cross-validation ----------------------------------------
keep_models <- c(intercept_only = "Intercept", activity = "Activity",
                 environment = "Environment", source = "Source",
                 all_domains_controls = "All domains")
MODEL_COLS <- c(Intercept = "#999999", Activity = "#56B4E9", Environment = "#E69F00",
                Source = "#0072B2", "All domains" = "#CC79A7")
cvd <- cv %>% filter(model %in% names(keep_models)) %>%
  mutate(model_lab = factor(keep_models[model], levels = keep_models),
         outcome_label = recode(outcome_label, "LAeq_24h" = "LAeq 24h"),
         outcome_label = factor(outcome_label,
                                levels = c("Linter-peak", "Lpeak", "Lnight", "Lday", "LAeq 24h")))
write_csv(cvd, file.path(DATA_DIR, "fig4_domain_cv.csv"))

pd <- ggplot(cvd, aes(cv_rmse_mean, outcome_label)) +
  geom_line(aes(group = outcome_label), colour = "grey80", linewidth = 0.4) +
  geom_point(aes(colour = model_lab, size = model_lab == "Source")) +
  scale_y_discrete(labels = c(
    expression(italic(L)[plain("inter-peak")]),
    expression(italic(L)[plain("peak")]),
    expression(italic(L)[plain("night")]),
    expression(italic(L)[plain("day")]),
    expression(italic(L)[plain("Aeq")]*" 24h"))) +
  scale_size_manual(values = c(`TRUE` = 2.2, `FALSE` = 1.5), guide = "none") +
  scale_colour_manual(values = MODEL_COLS, name = NULL) +
  scale_x_continuous(breaks = seq(3.5, 5.5, 0.5)) +
  labs(x = "Cross-validated RMSE (dB)", y = NULL) +
  theme_pub(BS, AT) +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.02),
        legend.justification = c(1, 0), legend.text = element_text(size = 6.8),
        legend.key.size = unit(0.30, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA))

# ---- Compose 2x2 (single grid: equal cells, aligned axes) --------------------
fig4 <- pa + pb + pc + pd +
  plot_layout(ncol = 2, widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY))

save_fig(fig4, "fig4_source_mechanism.png", width_mode = "double_column", height_in = 5.2)
message("Fig 4 done. n_veh=", n_veh,
        " rho_traffic=", round(get_rho("laeq_24h_leq_db_a","sound_traffic_share"),2),
        " rho_natural=", round(get_rho("laeq_24h_leq_db_a","sound_natural_share"),2))
