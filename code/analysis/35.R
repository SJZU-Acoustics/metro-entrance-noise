suppressPackageStartupMessages(library(tidyverse))
script <- "code/analysis/35.R"
analysis_root <- file.path("intermediate", "analysis_35_assessment_information")
exploration_root <- "intermediate"
source(file.path("code", "shared", "helpers.R"))
d <- load_p11(exploration_root)
out <- file.path(analysis_root, "outcome")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv(d$input_hashes, file.path(out, "input_sha256.csv"))
capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))

n <- nrow(d$master)
z <- scale(d$env)
hours <- 0:23
day <- hours %in% 7:22
night <- !day
L24 <- d$leq(d$levels)
Lday <- d$leq(d$levels[, day])
Lnight <- d$leq(d$levels[, night])
stopifnot(max(abs(Lday - d$noise$lday_leq_db_a)) < 1e-8,
          max(abs(Lnight - d$noise$lnight_leq_db_a)) < 1e-8)

st <- tibble(station_key = d$master$station_key, station_name = d$master$station_name_cn,
  date = d$master$measurement_date, month = format(as.Date(d$master$measurement_date), "%m"),
  day_type = d$master$day_type, L24 = L24, Lday = Lday, Lnight = Lnight,
  day_night_contrast = Lday - Lnight,
  traffic_hours = rowSums(d$groups == "traffic"), natural_hours = rowSums(d$groups == "natural"),
  human_hours = rowSums(d$groups == "human"), mechanical_hours = rowSums(d$groups == "mechanical"),
  natural_day_hours = rowSums(d$groups[, day] == "natural"),
  natural_night_hours = rowSums(d$groups[, night] == "natural"),
  daily_winner = d$sources[max.col(d$source, ties.method = "first")],
  winner_tied = apply(d$source, 1, function(x) sum(x == max(x)) > 1))
write_csv(st, file.path(out, "station_descriptors.csv"))

ij <- t(combn(n, 2)); i <- ij[, 1]; j <- ij[, 2]
env_diff <- abs(z[i, ] - z[j, ])
tv <- rowSums(abs(d$source[i, ] - d$source[j, ])) / 2
pairs <- tibble(i = i, j = j, a = st$station_key[i], b = st$station_key[j],
  level_gap = abs(L24[i] - L24[j]),
  day_gap = abs(Lday[i] - Lday[j]), night_gap = abs(Lnight[i] - Lnight[j]),
  period_max_gap = pmax(day_gap, night_gap),
  context_rms_sd = sqrt(rowMeans(env_diff^2)), context_max_sd = apply(env_diff, 1, max),
  same_month = st$month[i] == st$month[j], same_date = st$date[i] == st$date[j],
  same_daytype = st$day_type[i] == st$day_type[j],
  tv = tv, equivalent_reallocated_hours = 24 * tv,
  different_label_hours = rowSums(d$groups[i, ] != d$groups[j, ]),
  hourly_level_rmse = sqrt(rowMeans((d$levels[i, ] - d$levels[j, ])^2)),
  natural_hour_gap = abs(st$natural_hours[i] - st$natural_hours[j]),
  human_hour_gap = abs(st$human_hours[i] - st$human_hours[j]),
  different_daily_winner = st$daily_winner[i] != st$daily_winner[j] &
    !st$winner_tied[i] & !st$winner_tied[j])
write_csv(pairs, file.path(out, "all_pairs.csv"))

# Primary caliper: 1 dB; context similarity requires EACH of eight variables within 1 SD.
# Both are diagnostic choices, not equivalence tests, legal limits or treatment thresholds.
# Preserved sensitivity grid: 0.5/1/2 dB, context 0.5/1/1.5 SD, source TV 3/6/9 hours.
summ <- function(p, label) {
  tibble(comparison = label, pairs = nrow(p), stations = length(unique(c(p$a, p$b))),
    median_level_gap = median(p$level_gap), median_context_max_sd = median(p$context_max_sd),
    median_tv_hours = median(p$equivalent_reallocated_hours),
    tv_hours_q75 = unname(quantile(p$equivalent_reallocated_hours, .75)),
    tv_ge3 = sum(p$equivalent_reallocated_hours >= 3 - 1e-8),
    tv_ge6 = sum(p$equivalent_reallocated_hours >= 6 - 1e-8),
    tv_ge9 = sum(p$equivalent_reallocated_hours >= 9 - 1e-8),
    median_different_label_hours = median(p$different_label_hours),
    median_hourly_rmse = median(p$hourly_level_rmse),
    night_gap_ge3 = sum(p$night_gap >= 3), different_daily_winner = sum(p$different_daily_winner))
}
groups <- list(all_pairs = pairs,
  L24_1dB = filter(pairs, level_gap <= 1),
  L24_context = filter(pairs, level_gap <= 1, context_max_sd <= 1),
  L24_context_month_daytype = filter(pairs, level_gap <= 1, context_max_sd <= 1, same_month, same_daytype),
  L24_context_same_date = filter(pairs, level_gap <= 1, context_max_sd <= 1, same_date),
  day_night_1dB = filter(pairs, level_gap <= 1, period_max_gap <= 1),
  day_night_context = filter(pairs, level_gap <= 1, period_max_gap <= 1, context_max_sd <= 1),
  day_night_context_month_daytype = filter(pairs, level_gap <= 1, period_max_gap <= 1,
    context_max_sd <= 1, same_month, same_daytype),
  day_night_context_same_date = filter(pairs, level_gap <= 1, period_max_gap <= 1,
    context_max_sd <= 1, same_date))
summary <- imap_dfr(groups, summ)
write_csv(summary, file.path(out, "pair_summary.csv"))

# Disjoint greedy matching gives each station at most one appearance.
# For context matches, minimise context RMS first, then the permitted level gap.
# No maximum-contrast pair selection is used.
disjoint <- function(p, by_context = TRUE) {
  if (by_context) p <- arrange(p, context_rms_sd, level_gap, a, b)
  else p <- arrange(p, level_gap, a, b)
  used <- character(); keep <- integer()
  for (k in seq_len(nrow(p))) {
    if (!(p$a[k] %in% used || p$b[k] %in% used)) {
      keep <- c(keep, k); used <- c(used, p$a[k], p$b[k])
    }
  }
  p[keep, ]
}
matched <- imap(groups[-1], function(p, nm) disjoint(p, grepl("context", nm)))
write_csv(imap_dfr(matched, summ), file.path(out, "disjoint_summary.csv"))
write_csv(imap_dfr(matched, ~ mutate(.x, comparison = .y)), file.path(out, "disjoint_pairs.csv"))
