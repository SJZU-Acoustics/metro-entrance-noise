suppressPackageStartupMessages(library(tidyverse))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- "intermediate"
source(file.path("code", "shared", "helpers.R"))
d <- load_p11(root)
out <- file.path(root, "analysis_41_priority_profile_reference/outcome")
periods <- p11_periods()[c("day", "night")]; n <- nrow(d$master)
# Entrance-by-source counts of primary hourly winners in each period.
counts <- lapply(periods, function(h) sapply(d$sources, function(g) rowSums(d$groups[, h + 1, drop = FALSE] == g)))
share <- function(ids, p, g) sum(counts[[p]][ids, g]) / (length(ids) * length(periods[[p]]))
network <- map_dfr(names(periods), function(p) tibble(period = p, source = d$sources,
  hours = colSums(counts[[p]]), share = hours / sum(hours)))
stopifnot(all(abs(tapply(network$share, network$period, sum) - 1) < 1e-12))
s <- p11_scores(d$panel) %>% filter(period %in% names(periods)) %>% arrange(period, station_key)
rules <- c(joint = "raw50_joint", level = "level", activity = "raw50_activity")
rank_of <- function(rule, p) rank(-filter(s, period == p)[[rules[[rule]]]], ties.method = "first")
ranks <- lapply(names(rules), function(r) lapply(names(periods), function(p) rank_of(r, p)) %>% set_names(names(periods))) %>% set_names(names(rules))
# Four set-by-period cells and the symmetric decomposition, at every shortlist size and rule.
profile <- map_dfr(names(rules), function(r) map_dfr(1:30, function(k) {
  D <- which(ranks[[r]]$day <= k); N <- which(ranks[[r]]$night <= k)
  map_dfr(d$sources, function(g) tibble(rule = r, budget = k, source = g, union_size = length(union(D, N)),
    day_day = share(D, "day", g), day_night = share(D, "night", g), night_day = share(N, "day", g), night_night = share(N, "night", g)))
})) %>% mutate(total_difference = night_night - day_day,
  temporal_component = .5 * ((day_night - day_day) + (night_night - night_day)),
  membership_component = .5 * ((night_day - day_day) + (night_night - day_night))) %>%
  left_join(network %>% select(period, source, share) %>% pivot_wider(names_from = period, values_from = share, names_prefix = "network_"), by = "source")
stopifnot(max(abs(profile$total_difference - profile$temporal_component - profile$membership_component)) < 1e-12)
a38 <- read_csv(file.path(root, "analysis_38_priority_sources/outcome/source_decomposition.csv"), show_col_types = FALSE) %>%
  filter(sample == "all109", label == "source_group")
chk <- inner_join(filter(profile, rule == "joint"), a38, by = c("budget", "source"), suffix = c("", "_a38"))
stopifnot(nrow(chk) == 12, max(abs(chk$day_day - chk$day_day_a38)) < 1e-12, max(abs(chk$temporal_component - chk$temporal_component_a38)) < 1e-12)
# Random sets of the same size: a descriptive reference for how far a selected set departs from the network.
set.seed(20261002); B <- 10000L
random_reference <- map_dfr(1:30, function(k) {
  draws <- replicate(B, sample.int(n, k))
  if (k == 1) draws <- matrix(draws, nrow = 1)
  map_dfr(names(periods), function(p) map_dfr(c("traffic", "human"), function(g) {
    v <- colSums(matrix(counts[[p]][draws, g], nrow = k)) / (k * length(periods[[p]]))
    obs <- share(which(ranks$joint[[p]] <= k), p, g)
    tibble(budget = k, period = p, source = g, observed_joint = obs, random_q025 = quantile(v, .025), random_median = median(v),
      random_q975 = quantile(v, .975), random_fraction_at_or_below = mean(v <= obs), random_fraction_at_or_above = mean(v >= obs))
  }))
})
# Own-period shares by rank band, for joint and level ranks.
bands <- tibble(band = c("1-5", "6-10", "11-20", "21-30", "31-109"), lo = c(1, 6, 11, 21, 31), hi = c(5, 10, 20, 30, 109))
rank_bands <- map_dfr(c("joint", "level"), function(r) map_dfr(names(periods), function(p) pmap_dfr(bands, function(band, lo, hi) {
  ids <- which(ranks[[r]][[p]] >= lo & ranks[[r]][[p]] <= hi)
  tibble(rule = r, period = p, band = band, entrances = length(ids), source = d$sources,
    share = map_dbl(d$sources, ~share(ids, p, .x)), sites_present = map_int(d$sources, ~sum(counts[[p]][ids, .x] > 0)))
})))
# Leave-one-entrance-out range of the own-period traffic share, without reselection.
loo <- map_dfr(c(5, 10, 20), function(k) map_dfr(names(periods), function(p) {
  ids <- which(ranks$joint[[p]] <= k); v <- map_dbl(ids, ~share(setdiff(ids, .x), p, "traffic"))
  tibble(budget = k, period = p, traffic_share = share(ids, p, "traffic"), loo_min = min(v), loo_max = max(v))
}))
p11_save(list(network_reference = network, profile_by_size = profile, random_reference = random_reference,
  rank_bands = rank_bands, traffic_leave_one_out = loo), out, d, script, root)
print(network, n = Inf)
print(profile %>% filter(rule == "joint", source == "traffic", budget %in% c(5, 10, 15, 20, 30)) %>%
  mutate(across(day_day:membership_component, ~round(100 * .x, 2))) %>% select(budget, union_size, day_day:membership_component), width = Inf)
print(profile %>% filter(source == "traffic", budget %in% c(5, 10, 20)) %>% transmute(rule, budget, DD = round(100 * day_day, 1), NN = round(100 * night_night, 1)) %>% arrange(budget, rule), n = Inf)
print(random_reference %>% filter(budget %in% c(5, 10, 15, 20, 30)) %>% mutate(across(observed_joint:random_q975, ~round(100 * .x, 1))), n = Inf, width = Inf)
print(rank_bands %>% filter(source %in% c("traffic", "human")) %>% mutate(share = round(100 * share, 1)) %>% pivot_wider(names_from = c(period, source), values_from = c(share, sites_present)), n = Inf, width = Inf)
print(loo, width = Inf)
