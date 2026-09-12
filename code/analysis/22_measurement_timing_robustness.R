suppressPackageStartupMessages({
  library(tidyverse)
})

# Paths are relative to the repository root; python run_all.py runs every module there.
data_dir <- file.path("intermediate", "clean")
out_dir <- file.path("intermediate", "analysis")
if (!dir.exists(data_dir)) stop("Run from the repository root after code/load_data.py has built intermediate/clean/", call. = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(0.10, "cm"),
      axis.title = element_text(size = axis_title_size, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      legend.text = element_text(size = base_size),
      legend.title = element_blank(),
      legend.key = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

extract_coef <- function(model, outcome, model_name) {
  sm <- summary(model)
  coef_tbl <- as.data.frame(sm$coefficients)
  coef_tbl$term <- rownames(coef_tbl)
  rownames(coef_tbl) <- NULL
  coef_tbl %>%
    as_tibble() %>%
    rename(coef = Estimate, se = `Std. Error`, t_value = `t value`, p_value = `Pr(>|t|)`) %>%
    filter(term %in% c("sound_traffic_share", "sound_natural_share", "sound_human_share")) %>%
    mutate(outcome = outcome, model = model_name) %>%
    select(outcome, model, term, coef, se, t_value, p_value)
}

noise <- read_csv(file.path(data_dir, "noise_metrics_station.csv"), show_col_types = FALSE) %>%
  select(station_key, laeq_24h_leq_db_a, lday_leq_db_a, lnight_leq_db_a, lpeak_leq_db_a, linter_peak_leq_db_a)

master <- read_csv(file.path(data_dir, "station_master.csv"), show_col_types = FALSE) %>%
  mutate(
    poi_total_density_300m = rowSums(across(ends_with("_poi_density_300m")), na.rm = TRUE),
    measurement_date = as.Date(measurement_date),
    measurement_month = format(measurement_date, "%Y-%m"),
    measurement_month_num = as.integer(format(measurement_date, "%m")),
    measurement_weekday = weekdays(measurement_date),
    line_name = factor(line_name),
    station_type_code = factor(station_type_code),
    day_type = factor(day_type),
    measurement_month = factor(measurement_month)
  )

data <- master %>%
  select(
    station_key, station_name_cn, line_name, station_type_code, station_type_name_cn,
    measurement_date, measurement_month, measurement_month_num, measurement_weekday, day_type,
    sound_traffic_share, sound_natural_share, sound_human_share,
    `100m_urban_arterial_road_density_m_m`,
    `500m_parks_green_space_land_pct`,
    `500m_residential_land_pct`,
    `500m_commercial_services_land_pct`,
    `300m_parks_green_space_land_pct`,
    `100m_shannon_entropy_index`,
    poi_total_density_300m,
    `500_ndvi`
  ) %>%
  left_join(noise, by = "station_key") %>%
  filter(!is.na(laeq_24h_leq_db_a))

outcomes <- c("laeq_24h_leq_db_a", "lday_leq_db_a", "lnight_leq_db_a", "lpeak_leq_db_a", "linter_peak_leq_db_a")
outcome_labels <- c(
  laeq_24h_leq_db_a = "LAeq_24h",
  lday_leq_db_a = "Lday",
  lnight_leq_db_a = "Lnight",
  lpeak_leq_db_a = "Lpeak",
  linter_peak_leq_db_a = "Linter-peak"
)

distribution <- data %>%
  count(measurement_month, day_type, name = "n_stations") %>%
  arrange(measurement_month, day_type)
write_csv(distribution, file.path(out_dir, "22_timing__sample_distribution_month_daytype.csv"))

line_month_distribution <- data %>%
  count(line_name, measurement_month, name = "n_stations") %>%
  arrange(line_name, measurement_month)
write_csv(line_month_distribution, file.path(out_dir, "22_timing__sample_distribution_line_month.csv"))

month_kw <- map_dfr(outcomes, function(outcome) {
  d <- data %>% select(measurement_month, value = all_of(outcome)) %>% drop_na()
  test <- kruskal.test(value ~ measurement_month, data = d)
  tibble(
    outcome = outcome,
    outcome_label = outcome_labels[[outcome]],
    n = nrow(d),
    k_months = n_distinct(d$measurement_month),
    statistic = unname(test$statistic),
    p_value = test$p.value
  )
})
write_csv(month_kw, file.path(out_dir, "22_timing__month_kruskal_wallis.csv"))

daytype_tests <- map_dfr(outcomes, function(outcome) {
  d <- data %>% select(day_type, value = all_of(outcome)) %>% drop_na()
  test <- wilcox.test(value ~ day_type, data = d, exact = FALSE)
  stats <- d %>%
    group_by(day_type) %>%
    summarise(n = n(), mean = mean(value), median = median(value), .groups = "drop")
  tibble(
    outcome = outcome,
    outcome_label = outcome_labels[[outcome]],
    n_weekday = stats$n[stats$day_type == "weekday"][1],
    n_weekend = stats$n[stats$day_type == "weekend"][1],
    median_weekday = stats$median[stats$day_type == "weekday"][1],
    median_weekend = stats$median[stats$day_type == "weekend"][1],
    median_weekend_minus_weekday = median_weekend - median_weekday,
    p_value = test$p.value
  )
})
write_csv(daytype_tests, file.path(out_dir, "22_timing__daytype_wilcox.csv"))

env_terms <- paste(c(
  "`100m_urban_arterial_road_density_m_m`",
  "`500m_parks_green_space_land_pct`",
  "`500m_residential_land_pct`",
  "`500m_commercial_services_land_pct`",
  "`300m_parks_green_space_land_pct`",
  "`100m_shannon_entropy_index`",
  "poi_total_density_300m",
  "`500_ndvi`"
), collapse = " + ")
source_terms <- "sound_traffic_share + sound_natural_share + sound_human_share"
base_rhs <- paste(source_terms, env_terms, "line_name + station_type_code + day_type", sep = " + ")
month_rhs <- paste(base_rhs, "measurement_month", sep = " + ")

coef_rows <- list()
increment_rows <- list()

for (outcome in outcomes) {
  base_formula <- as.formula(paste(outcome, "~", base_rhs))
  month_formula <- as.formula(paste(outcome, "~", month_rhs))
  model_base <- lm(base_formula, data = data)
  model_month <- lm(month_formula, data = data)
  coef_rows[[length(coef_rows) + 1]] <- extract_coef(model_base, outcome, "base_controls")
  coef_rows[[length(coef_rows) + 1]] <- extract_coef(model_month, outcome, "month_added")
  comp <- anova(model_base, model_month)
  increment_rows[[length(increment_rows) + 1]] <- tibble(
    outcome = outcome,
    outcome_label = outcome_labels[[outcome]],
    n = nobs(model_month),
    base_adj_r2 = summary(model_base)$adj.r.squared,
    month_adj_r2 = summary(model_month)$adj.r.squared,
    delta_adj_r2 = month_adj_r2 - base_adj_r2,
    month_increment_p = comp$`Pr(>F)`[2]
  )
}

coef_df <- bind_rows(coef_rows) %>%
  group_by(outcome, term) %>%
  mutate(
    base_coef = coef[model == "base_controls"][1],
    month_coef = coef[model == "month_added"][1],
    sign_changed = sign(base_coef) != sign(month_coef)
  ) %>%
  ungroup()
increment_df <- bind_rows(increment_rows)

write_csv(coef_df, file.path(out_dir, "22_timing__source_coefficients_with_month.csv"))
write_csv(increment_df, file.path(out_dir, "22_timing__month_incremental_tests.csv"))

summary_tbl <- tibble(
  item = c(
    "n_valid_stations",
    "n_months",
    "n_weekday",
    "n_weekend",
    "month_kw_sig_0_05",
    "daytype_sig_0_05",
    "source_sign_changes_after_month",
    "month_increment_sig_0_05"
  ),
  value = c(
    nrow(data),
    n_distinct(data$measurement_month),
    sum(data$day_type == "weekday"),
    sum(data$day_type == "weekend"),
    sum(month_kw$p_value < 0.05),
    sum(daytype_tests$p_value < 0.05),
    sum(coef_df %>% filter(model == "month_added") %>% pull(sign_changed), na.rm = TRUE),
    sum(increment_df$month_increment_p < 0.05, na.rm = TRUE)
  )
)
write_csv(summary_tbl, file.path(out_dir, "22_timing__summary.csv"))

p_dist <- distribution %>%
  ggplot(aes(measurement_month, n_stations, fill = day_type)) +
  geom_col(position = "stack", width = 0.72) +
  scale_fill_manual(values = c(weekday = "#0072B2", weekend = "#E69F00")) +
  labs(x = "Measurement month", y = "Stations (n)") +
  theme_pub() +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "22_timing__sample_distribution_month_daytype.png"), p_dist, width = 85 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

p_coef <- coef_df %>%
  mutate(
    term = recode(term, sound_traffic_share = "Traffic", sound_natural_share = "Natural", sound_human_share = "Human"),
    outcome_label = factor(outcome_labels[outcome], levels = outcome_labels[outcomes]),
    model = recode(model, base_controls = "Base", month_added = "Month added")
  ) %>%
  ggplot(aes(model, coef, group = term, colour = term)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.2) +
  facet_wrap(~outcome_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(Traffic = "#D55E00", Natural = "#009E73", Human = "#0072B2")) +
  labs(x = NULL, y = "Coefficient") +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "22_timing__source_coefficients_with_month.png"), p_coef, width = 178 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

message("Completed measurement timing robustness for ", nrow(data), " stations.")
