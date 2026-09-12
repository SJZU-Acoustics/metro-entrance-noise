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

period_hours <- list(
  day = 7:22,
  night = c(23, 0:6),
  peak = c(7, 8, 17, 18),
  offpeak = 10:16
)

safe_log1p <- function(x) log1p(pmax(as.numeric(x), 0))

make_folds <- function(n, k = 5, repeats = 100, seed = 20260529) {
  set.seed(seed)
  folds <- vector("list", repeats * k)
  idx <- 1
  for (r in seq_len(repeats)) {
    fold_id <- sample(rep(seq_len(k), length.out = n))
    for (fold in seq_len(k)) {
      folds[[idx]] <- list(rep_id = r, fold = fold, test = which(fold_id == fold))
      idx <- idx + 1
    }
  }
  folds
}

fit_predict_lm <- function(x, y, train_idx, test_idx) {
  x_train <- x[train_idx, , drop = FALSE]
  x_test <- x[test_idx, , drop = FALSE]
  keep <- apply(x_train, 2, function(col) sd(col, na.rm = TRUE) > 0)
  x_train <- x_train[, keep, drop = FALSE]
  x_test <- x_test[, keep, drop = FALSE]
  design_train <- cbind("(Intercept)" = 1, x_train)
  design_test <- cbind("(Intercept)" = 1, x_test)
  fit <- lm.fit(design_train, y[train_idx])
  coef <- fit$coefficients
  coef[is.na(coef)] <- 0
  as.numeric(design_test %*% coef)
}

cv_rmse <- function(x, y, folds) {
  preds <- rep(NA_real_, length(y))
  records <- vector("list", length(folds))
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]$test
    train_idx <- setdiff(seq_along(y), test_idx)
    pred <- fit_predict_lm(x, y, train_idx, test_idx)
    records[[i]] <- tibble(
      rep_id = folds[[i]]$rep_id,
      fold = folds[[i]]$fold,
      rmse = sqrt(mean((y[test_idx] - pred)^2, na.rm = TRUE))
    )
    preds[test_idx] <- pred
  }
  bind_rows(records)
}

full_fit_stats <- function(x, y) {
  keep <- apply(x, 2, function(col) sd(col, na.rm = TRUE) > 0)
  x <- x[, keep, drop = FALSE]
  design <- cbind("(Intercept)" = 1, x)
  fit <- lm.fit(design, y)
  fitted <- as.numeric(design %*% replace(fit$coefficients, is.na(fit$coefficients), 0))
  rss <- sum((y - fitted)^2, na.rm = TRUE)
  tss <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- 1 - rss / tss
  p <- fit$rank - 1
  n <- length(y)
  adj_r2 <- 1 - (1 - r2) * (n - 1) / max(n - p - 1, 1)
  tibble(n = n, p = p, r2 = r2, adj_r2 = adj_r2, full_rmse = sqrt(mean((y - fitted)^2, na.rm = TRUE)))
}

noise <- read_csv(file.path(data_dir, "noise_metrics_station.csv"), show_col_types = FALSE) %>%
  select(station_key, laeq_24h_leq_db_a, lday_leq_db_a, lnight_leq_db_a, lpeak_leq_db_a, linter_peak_leq_db_a)
master <- read_csv(file.path(data_dir, "station_master.csv"), show_col_types = FALSE) %>%
  mutate(poi_total_density_300m = rowSums(across(ends_with("_poi_density_300m")), na.rm = TRUE))
panel <- read_csv(file.path(data_dir, "station_hour_panel.csv"), show_col_types = FALSE) %>%
  mutate(pop_heat_50m = as.numeric(pop_heat_50m))

heat_features <- imap_dfr(period_hours, function(hours, period_name) {
  panel %>%
    filter(hour %in% hours) %>%
    group_by(station_key) %>%
    summarise(value = mean(pop_heat_50m, na.rm = TRUE), .groups = "drop") %>%
    mutate(period_feature = paste0("heat_", period_name, "_mean"))
}) %>%
  pivot_wider(names_from = period_feature, values_from = value)

data <- master %>%
  select(
    station_key, line_name, station_type_code, day_type,
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
  left_join(heat_features, by = "station_key") %>%
  left_join(noise, by = "station_key") %>%
  mutate(
    across(starts_with("heat_"), safe_log1p),
    line_name = factor(line_name),
    station_type_code = factor(station_type_code),
    day_type = factor(day_type)
  )

outcomes <- c("laeq_24h_leq_db_a", "lday_leq_db_a", "lnight_leq_db_a", "lpeak_leq_db_a", "linter_peak_leq_db_a")
outcome_labels <- c(
  laeq_24h_leq_db_a = "LAeq_24h",
  lday_leq_db_a = "Lday",
  lnight_leq_db_a = "Lnight",
  lpeak_leq_db_a = "Lpeak",
  linter_peak_leq_db_a = "Linter-peak"
)
env_vars <- c(
  "`100m_urban_arterial_road_density_m_m`",
  "`500m_parks_green_space_land_pct`",
  "`500m_residential_land_pct`",
  "`500m_commercial_services_land_pct`",
  "`300m_parks_green_space_land_pct`",
  "`100m_shannon_entropy_index`",
  "poi_total_density_300m",
  "`500_ndvi`"
)
source_vars <- c("sound_traffic_share", "sound_natural_share", "sound_human_share")
control_vars <- c("line_name", "station_type_code", "day_type")
activity_map <- c(
  laeq_24h_leq_db_a = "heat_day_mean + heat_night_mean",
  lday_leq_db_a = "heat_day_mean",
  lnight_leq_db_a = "heat_night_mean",
  lpeak_leq_db_a = "heat_peak_mean",
  linter_peak_leq_db_a = "heat_offpeak_mean"
)

model_terms <- function(outcome) {
  activity <- activity_map[[outcome]]
  list(
    intercept_only = "1",
    controls = paste(control_vars, collapse = " + "),
    environment = paste(env_vars, collapse = " + "),
    source = paste(source_vars, collapse = " + "),
    activity = activity,
    environment_source = paste(c(env_vars, source_vars), collapse = " + "),
    exposure_domains = paste(c(env_vars, source_vars, activity), collapse = " + "),
    all_domains_controls = paste(c(control_vars, env_vars, source_vars, activity), collapse = " + ")
  )
}

performance <- list()
fold_records <- list()

for (outcome in outcomes) {
  terms <- model_terms(outcome)
  needed <- unique(unlist(str_extract_all(paste(terms, collapse = " + "), "`[^`]+`|[A-Za-z0-9_]+")))
  needed <- str_remove_all(needed, "`")
  needed <- setdiff(needed, c("1"))
  model_data <- data %>%
    select(all_of(c(outcome, needed))) %>%
    drop_na()
  y <- model_data[[outcome]]
  folds <- make_folds(length(y), k = 5, repeats = 100)
  for (model_name in names(terms)) {
    rhs <- terms[[model_name]]
    if (rhs == "1") {
      x <- matrix(numeric(0), nrow = length(y), ncol = 0)
    } else {
      x <- model.matrix(as.formula(paste("~", rhs)), data = model_data)[, -1, drop = FALSE]
    }
    cv <- cv_rmse(x, y, folds)
    stats <- full_fit_stats(x, y)
    performance[[length(performance) + 1]] <- stats %>%
      mutate(
        outcome = outcome,
        outcome_label = outcome_labels[[outcome]],
        model = model_name,
        cv_rmse_mean = mean(cv$rmse),
        cv_rmse_sd = sd(cv$rmse),
        cv_rmse_median = median(cv$rmse)
      )
    fold_records[[length(fold_records) + 1]] <- cv %>%
      mutate(outcome = outcome, model = model_name)
  }
}

perf <- bind_rows(performance) %>%
  relocate(outcome, outcome_label, model)
folds_out <- bind_rows(fold_records)
write_csv(perf, file.path(out_dir, "21_domain_cv__model_performance.csv"))
write_csv(folds_out, file.path(out_dir, "21_domain_cv__fold_rmse.csv"))

# Paired comparison across the 100 repeats: per repeat, the mean RMSE over the
# five folds for each model; count the repeats in which the source-only model
# has the lower error than each rival model.
rep_means <- folds_out %>%
  group_by(outcome, model, rep_id) %>%
  summarise(rmse = mean(rmse), .groups = "drop")
paired <- rep_means %>%
  filter(model != "source") %>%
  left_join(rep_means %>% filter(model == "source") %>% select(outcome, rep_id, source_rmse = rmse),
            by = c("outcome", "rep_id")) %>%
  group_by(outcome, model) %>%
  summarise(
    n_repeats = n(),
    source_wins = sum(source_rmse < rmse),
    mean_diff_vs_source = mean(rmse - source_rmse),
    .groups = "drop"
  ) %>%
  mutate(outcome_label = outcome_labels[outcome]) %>%
  relocate(outcome, outcome_label, model)
write_csv(paired, file.path(out_dir, "21_domain_cv__paired_wins.csv"))

gain <- perf %>%
  select(outcome, outcome_label, model, cv_rmse_mean, adj_r2) %>%
  group_by(outcome) %>%
  mutate(
    null_cv_rmse = cv_rmse_mean[model == "intercept_only"][1],
    controls_cv_rmse = cv_rmse_mean[model == "controls"][1],
    cv_rmse_gain_vs_null = null_cv_rmse - cv_rmse_mean,
    cv_rmse_gain_vs_controls = controls_cv_rmse - cv_rmse_mean,
    best_model = model[which.min(cv_rmse_mean)][1],
    rank_cv = rank(cv_rmse_mean, ties.method = "first")
  ) %>%
  ungroup()
write_csv(gain, file.path(out_dir, "21_domain_cv__model_gain.csv"))

summary_tbl <- gain %>%
  filter(rank_cv == 1) %>%
  count(best_model, name = "n_outcomes_best") %>%
  arrange(desc(n_outcomes_best))
write_csv(summary_tbl, file.path(out_dir, "21_domain_cv__summary.csv"))

model_levels <- c("intercept_only", "controls", "activity", "environment", "source", "environment_source", "exposure_domains", "all_domains_controls")
model_labels <- c(
  intercept_only = "Intercept only",
  controls = "Controls",
  activity = "Activity",
  environment = "Environment",
  source = "Sources",
  environment_source = "Env + sources",
  exposure_domains = "Activity + env + sources",
  all_domains_controls = "All + controls"
)

p_perf <- perf %>%
  mutate(
    model = factor(model, levels = model_levels),
    model_label = recode(as.character(model), !!!model_labels),
    model_label = factor(model_label, levels = model_labels),
    outcome_label = factor(outcome_label, levels = outcome_labels[outcomes])
  ) %>%
  ggplot(aes(model_label, cv_rmse_mean, colour = outcome_label, group = outcome_label)) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7")) +
  labs(x = NULL, y = "Repeated 5-fold CV RMSE (dB)") +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top"
  )
ggsave(file.path(out_dir, "21_domain_cv__model_performance.png"), p_perf, width = 178 / 25.4, height = 105 / 25.4, dpi = 600, bg = "white")

message("Completed domain CV models.")
