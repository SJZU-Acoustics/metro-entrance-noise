suppressPackageStartupMessages(library(tidyverse))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- "intermediate"
source(file.path("code", "shared", "helpers.R"))
d <- load_p11(root)
out <- file.path(root, "analysis_39_persistent_error/outcome")
periods <- p11_periods()[c("day", "peak", "offpeak", "night")]
n <- nrow(d$master); B <- 1000L; set.seed(20261001)
A <- matrix(d$panel$pop_heat_50m, nrow = n, byrow = TRUE)
L <- d$levels; original <- p11_scores(d$panel)
scenarios <- crossing(scale = c("main", "stress"), persistence = c(0, .5, 1)) %>%
  mutate(level_sd = if_else(scale == "main", 1, 2), activity_cv = if_else(scale == "main", .2, .3),
    scenario = paste0(scale, "_rho", persistence))
replicates <- list(); frequencies <- list()
for (q in seq_len(nrow(scenarios))) {
  cfg <- scenarios[q, ]; counts <- matrix(0, nrow = n, ncol = length(periods)); rr <- list(); ix <- 0L
  for (b in seq_len(B)) {
    zL <- sqrt(cfg$persistence) * matrix(rnorm(n), n, 24) + sqrt(1 - cfg$persistence) * matrix(rnorm(n * 24), n, 24)
    zA <- sqrt(cfg$persistence) * matrix(rnorm(n), n, 24) + sqrt(1 - cfg$persistence) * matrix(rnorm(n * 24), n, 24)
    log_sd <- sqrt(log(1 + cfg$activity_cv^2))
    E <- 10^((L + cfg$level_sd * zL) / 10)
    Ap <- A * exp(log_sd * zA - log_sd^2 / 2)
    dayv <- rowSums((Ap * E)[, periods$day + 1, drop = FALSE])
    day_sel <- order(-dayv, seq_len(n))[1:10]
    pp <- 0L
    for (p in names(periods)) {
      pp <- pp + 1L; h <- periods[[p]] + 1
      joint <- rowSums((Ap * E)[, h, drop = FALSE]); product <- rowSums(Ap[, h, drop = FALSE]) * rowMeans(E[, h, drop = FALSE])
      lev <- rowMeans(E[, h, drop = FALSE])
      select <- order(-joint, seq_len(n))[1:10]; prod_sel <- order(-product, seq_len(n))[1:10]; lev_sel <- order(-lev, seq_len(n))[1:10]
      obs <- filter(original, period == p) %>% arrange(station_key)
      stopifnot(identical(obs$station_key, d$master$station_key))
      obs_sel <- order(-obs$raw50_joint, seq_len(n))[1:10]
      counts[select, pp] <- counts[select, pp] + 1
      ix <- ix + 1L
      rr[[ix]] <- tibble(scenario = cfg$scenario, replicate = b, period = p,
        observed_overlap = sum(select %in% obs_sel), product_overlap = sum(select %in% prod_sel),
        day_overlap = sum(select %in% day_sel),
        product_gap_pp = 100 * (sum(joint[select]) - sum(joint[prod_sel])) / sum(joint),
        day_gap_pp = 100 * (sum(joint[select]) - sum(joint[day_sel])) / sum(joint),
        level_gap_pp = 100 * (sum(joint[select]) - sum(joint[lev_sel])) / sum(joint),
        observed_list_gap_pp = 100 * (sum(joint[select]) - sum(joint[obs_sel])) / sum(joint))
    }
  }
  replicates[[q]] <- bind_rows(rr)
  frequencies[[q]] <- map_dfr(seq_along(periods), function(i) tibble(scenario = cfg$scenario,
    period = names(periods)[i], station_key = d$master$station_key, inclusion_frequency = counts[, i] / B))
}
reps <- bind_rows(replicates); freq <- bind_rows(frequencies)
summary <- reps %>% group_by(scenario, period) %>% summarise(across(c(observed_overlap, product_overlap,
  product_gap_pp, day_gap_pp, level_gap_pp, observed_list_gap_pp), list(mean = mean, q05 = ~quantile(.x, .05), q95 = ~quantile(.x, .95)),
  .names = "{.col}_{.fn}"), .groups = "drop") %>% left_join(scenarios, by = "scenario")
p11_save(list(scenarios = scenarios, simulation_summary = summary, inclusion_frequencies = freq,
  simulation_replicates = reps), out, d, script, root)
print(summary %>% filter(period == "night") %>% select(scenario, observed_overlap_mean, observed_overlap_q05,
  product_overlap_mean, product_gap_pp_mean, day_gap_pp_mean, observed_list_gap_pp_mean), width = Inf)
