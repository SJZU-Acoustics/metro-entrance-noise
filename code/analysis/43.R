suppressPackageStartupMessages(library(tidyverse))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- "intermediate"
source(file.path("code", "shared", "helpers.R"))
d <- load_p11(root)
out <- file.path(root, "analysis_43_recording_dates/outcome")
n <- nrow(d$master); periods <- p11_periods()[c("day", "night")]
dates <- d$master %>% count(measurement_date, name = "entrances") %>% mutate(date_id = row_number())
date_id <- dates$date_id[match(d$master$measurement_date, dates$measurement_date)]
stopifnot(nrow(dates) == 25, sum(dates$entrances) == n, !anyNA(date_id))
A <- matrix(d$panel$pop_heat_50m, nrow = n, byrow = TRUE); L <- d$levels; E <- 10^(L / 10)
rk <- function(v) rank(-v, ties.method = "first")
traffic <- lapply(periods, function(h) rowMeans(d$groups[, h + 1, drop = FALSE] == "traffic"))
# Headline selection and source quantities for any subset of entrances, reselecting within it.
headline <- function(keep) {
  J <- lapply(c(periods, list(full = 0:23)), function(h) rowSums((A * E)[keep, h + 1, drop = FALSE]))
  P <- lapply(periods, function(h) rowSums(A[keep, h + 1, drop = FALSE]) * rowMeans(E[keep, h + 1, drop = FALSE]))
  lev <- lapply(periods, function(h) rowMeans(E[keep, h + 1, drop = FALSE]))
  r <- lapply(J, rk); cv <- function(sel, p) 100 * sum(J[[p]][sel]) / sum(J[[p]])
  tibble(entrances = length(keep),
    night_gain_over_full = cv(r$night <= 10, "night") - cv(r$full <= 10, "night"),
    night_shared_with_full = sum(r$night <= 10 & r$full <= 10), night_shared_with_day = sum(r$night <= 10 & r$day <= 10),
    night_shared_with_level = sum(r$night <= 10 & rk(lev$night) <= 10), day_shared_with_level = sum(r$day <= 10 & rk(lev$day) <= 10),
    night_product_gap = cv(r$night <= 10, "night") - cv(rk(P$night) <= 10, "night"),
    day_product_gap = cv(r$day <= 10, "day") - cv(rk(P$day) <= 10, "day"),
    day_traffic_share = 100 * mean(traffic$day[keep][r$day <= 10]), night_traffic_share = 100 * mean(traffic$night[keep][r$night <= 10]),
    day_list = paste(sort(d$master$station_key[keep][r$day <= 10]), collapse = " "),
    night_list = paste(sort(d$master$station_key[keep][r$night <= 10]), collapse = " "))
}
full_sample <- headline(seq_len(n))
orig_day <- strsplit(full_sample$day_list, " ")[[1]]; orig_night <- strsplit(full_sample$night_list, " ")[[1]]
omit_date <- map_dfr(dates$date_id, function(k) {
  keep <- which(date_id != k); h <- headline(keep)
  removed <- d$master$station_key[date_id == k]
  h %>% mutate(measurement_date = dates$measurement_date[k], omitted_entrances = length(removed),
    omitted_day_priorities = sum(removed %in% orig_day), omitted_night_priorities = sum(removed %in% orig_night),
    day_retained = sum(strsplit(day_list, " ")[[1]] %in% orig_day), night_retained = sum(strsplit(night_list, " ")[[1]] %in% orig_night), .before = 1)
})
# Omit each entrance on any ten-entrance list (day, night or 24-hour) and reselect, with the same headline quantities.
full_rank <- rk(rowSums(A * E))
listed <- which(d$master$station_key %in% c(orig_day, orig_night) | full_rank <= 10)
omit_entrance <- map_dfr(listed, function(i) headline(setdiff(seq_len(n), i)) %>%
  mutate(omitted = d$master$station_key[i], day_retained = sum(strsplit(day_list, " ")[[1]] %in% orig_day),
    night_retained = sum(strsplit(night_list, " ")[[1]] %in% orig_night), .before = 1))
priority_dates <- d$master %>% transmute(station_key, measurement_date, day_priority = station_key %in% orig_day,
  night_priority = station_key %in% orig_night) %>% group_by(measurement_date) %>%
  summarise(entrances = n(), day_priorities = sum(day_priority), night_priorities = sum(night_priority), .groups = "drop")
# Persistent error shared by every entrance recorded on the same date, with the remainder independent by hour.
B <- 1000L; set.seed(20261002)
scenarios <- crossing(scale = c("main", "stress"), date_share = c(.5, 1)) %>%
  mutate(level_sd = if_else(scale == "main", 1, 2), activity_cv = if_else(scale == "main", .2, .3), scenario = paste0(scale, "_date", date_share))
obs_sel <- lapply(periods, function(h) order(-rowSums((A * E)[, h + 1, drop = FALSE]), seq_len(n))[1:10])
reps <- list(); freq <- list()
for (q in seq_len(nrow(scenarios))) {
  cfg <- scenarios[q, ]; counts <- matrix(0, n, length(periods)); rr <- vector("list", B * length(periods)); ix <- 0L
  log_sd <- sqrt(log(1 + cfg$activity_cv^2))
  for (b in seq_len(B)) {
    zL <- sqrt(cfg$date_share) * matrix(rnorm(25)[date_id], n, 24) + sqrt(1 - cfg$date_share) * matrix(rnorm(n * 24), n, 24)
    zA <- sqrt(cfg$date_share) * matrix(rnorm(25)[date_id], n, 24) + sqrt(1 - cfg$date_share) * matrix(rnorm(n * 24), n, 24)
    Ep <- 10^((L + cfg$level_sd * zL) / 10); Ap <- A * exp(log_sd * zA - log_sd^2 / 2)
    for (pp in seq_along(periods)) {
      h <- periods[[pp]] + 1; joint <- rowSums((Ap * Ep)[, h, drop = FALSE])
      product <- rowSums(Ap[, h, drop = FALSE]) * rowMeans(Ep[, h, drop = FALSE])
      sel <- order(-joint, seq_len(n))[1:10]; prod_sel <- order(-product, seq_len(n))[1:10]
      counts[sel, pp] <- counts[sel, pp] + 1; ix <- ix + 1L
      rr[[ix]] <- tibble(scenario = cfg$scenario, replicate = b, period = names(periods)[pp],
        observed_overlap = sum(sel %in% obs_sel[[pp]]),
        product_gap_pp = 100 * (sum(joint[sel]) - sum(joint[prod_sel])) / sum(joint),
        observed_list_gap_pp = 100 * (sum(joint[sel]) - sum(joint[obs_sel[[pp]]])) / sum(joint))
    }
  }
  reps[[q]] <- bind_rows(rr)
  freq[[q]] <- map_dfr(seq_along(periods), ~tibble(scenario = cfg$scenario, period = names(periods)[.x],
    station_key = d$master$station_key, inclusion_frequency = counts[, .x] / B))
}
reps <- bind_rows(reps); freq <- bind_rows(freq)
date_error_summary <- reps %>% group_by(scenario, period) %>% summarise(across(c(observed_overlap, product_gap_pp, observed_list_gap_pp),
  list(mean = mean, q05 = ~quantile(.x, .05), q95 = ~quantile(.x, .95)), .names = "{.col}_{.fn}"), .groups = "drop") %>%
  left_join(scenarios, by = "scenario")
p11_save(list(recording_dates = priority_dates, full_sample = full_sample, omit_one_date = omit_date, omit_one_entrance = omit_entrance,
  date_error_scenarios = scenarios, date_error_summary = date_error_summary, date_error_replicates = reps,
  date_error_inclusion = freq), out, d, script, root)
print(filter(priority_dates, day_priorities > 0 | night_priorities > 0), n = Inf)
print(select(full_sample, -day_list, -night_list), width = Inf)
print(omit_date %>% select(measurement_date, omitted_entrances, omitted_day_priorities, omitted_night_priorities, day_retained, night_retained,
  night_gain_over_full, night_shared_with_full, night_product_gap, day_traffic_share, night_traffic_share) %>% mutate(across(where(is.double), ~round(.x, 2))), n = Inf, width = Inf)
print(date_error_summary %>% select(scenario, period, observed_overlap_mean, observed_overlap_q05, observed_overlap_q95, product_gap_pp_mean, observed_list_gap_pp_mean), width = Inf)
