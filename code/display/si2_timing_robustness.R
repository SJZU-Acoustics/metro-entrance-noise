# SI Fig S4 -- Measurement-timing robustness of the source associations.
source(file.path("code", "display", "style.R"))

samp <- read_csv(file.path(EXPLORE_OUT, "22_timing__sample_distribution_month_daytype.csv"),
                 show_col_types = FALSE)
coef <- read_csv(file.path(EXPLORE_OUT, "22_timing__source_coefficients_with_month.csv"),
                 show_col_types = FALSE)
BS <- 8; AT <- 9

# Panel a: monitoring sample by month and day type
sa <- samp %>%
  mutate(month = factor(measurement_month),
         day_type = factor(str_to_title(day_type), levels = c("Weekday","Weekend")))
pa <- ggplot(sa, aes(month, n_stations, fill = day_type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(Weekday = "#0072B2", Weekend = "#E69F00"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Monitoring month", y = "Entrances") +
  theme_pub(BS, AT) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1), legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))

# Panel b: natural-share coefficient, base vs month-added, by outcome
out_lab <- c(laeq_24h_leq_db_a = "LAeq 24 h", lday_leq_db_a = "Lday",
             lnight_leq_db_a = "Lnight", lpeak_leq_db_a = "Lpeak",
             linter_peak_leq_db_a = "Linter-peak")
cb <- coef %>%
  filter(term == "sound_natural_share") %>%
  transmute(outcome = factor(out_lab[outcome], levels = rev(unname(out_lab))),
            model = recode(model, base_controls = "Base controls",
                           month_added = "+ month"),
            coef, lo = coef - 1.96 * se, hi = coef + 1.96 * se) %>%
  mutate(model = factor(model, levels = c("Base controls", "+ month")))

pb <- ggplot(cb, aes(coef, outcome, colour = model)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey55", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(size = 1.4, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Base controls" = "#009E73", "+ month" = "#CC79A7"),
                      name = NULL) +
  labs(x = "Natural-share coefficient (dB)", y = NULL) +
  theme_pub(BS, AT) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.02),
        legend.justification = c(0, 0), legend.text = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"))

fig <- pa + pb + plot_layout(widths = c(1, 1.1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY))
save_fig(fig, "si_timing_robustness.png", width_mode = "double_column", height_in = 2.9, dir = SIFIG_DIR)
message("SI S4 done.")
