# Independent arithmetic checks, including zero-error and common-shift invariants.
suppressPackageStartupMessages(library(tidyverse))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- "intermediate"
source("code/shared/helpers.R")
d <- load_p11(root)
readout <- function(a, f) read_csv(file.path(root, a, "outcome", paste0(f, ".csv")), show_col_types = FALSE)
a37 <- "analysis_37_information_choices"; a38 <- "analysis_38_priority_sources"; a39 <- "analysis_39_persistent_error"
cov <- readout(a37, "coverage"); ranks <- readout(a37, "rankings")
E <- 10^(d$levels / 10)
A <- matrix(d$panel$pop_heat_50m, nrow = 109, byrow = TRUE)
A20 <- matrix(d$panel$pop_heat_20m, nrow = 109, byrow = TRUE)
acts <- list(raw50 = A, raw20 = A20, sqrt50 = sqrt(A), cap50 = pmin(A, as.numeric(quantile(A, .95))))
periods <- p11_periods()
max_score_error <- 0; max_coverage_error <- 0
for (v in names(acts)) for (p in c("day", "peak", "offpeak", "night")) {
  h <- periods[[p]] + 1; hd <- periods$day + 1
  target <- rowSums((E * acts[[v]])[, h, drop = FALSE])
  baseline <- rowSums((E * A)[, h, drop = FALSE])
  vals <- list(full_level = rowMeans(E), period_level = rowMeans(E[, h, drop = FALSE]),
    full_joint = rowSums(E * acts[[v]]), period_joint = target,
    period_product = rowSums(acts[[v]][, h, drop = FALSE]) * rowMeans(E[, h, drop = FALSE]),
    day_joint = rowSums((E * acts[[v]])[, hd, drop = FALSE]),
    period_activity = rowSums(acts[[v]][, h, drop = FALSE]))
  for (m in names(vals)) {
    r <- ranks %>% filter(definition == v, period == p, method == m) %>% arrange(station_key)
    stopifnot(nrow(r) == 109, identical(r$station_key, d$master$station_key),
      identical(as.integer(r$rank), as.integer(rank(-vals[[m]], ties.method = "first"))))
    max_score_error <- max(max_score_error, max(abs(r$selection_score / vals[[m]] - 1)))
    z <- cov %>% filter(definition == v, period == p, method == m)
    for (i in seq_len(nrow(z))) {
      take <- r$rank <= z$budget[i]
      err <- max(abs(c(sum(target[take]) / sum(target) - z$fraction_index[i],
        sum(baseline[take]) / sum(baseline) - z$fraction_default_index[i])))
      max_coverage_error <- max(max_coverage_error, err)
    }
  }
}
stopifnot(nrow(cov) == 3360, nrow(ranks) == 12208, max_score_error < 1e-10, max_coverage_error < 1e-10)
src <- readout(a38, "station_profiles")
totals <- src %>% group_by(sample, station_key, period, label) %>% summarise(s = sum(share), h = sum(hours), .groups = "drop")
stopifnot(all(abs(totals$s - 1) < 1e-12), all(totals$h == lengths(periods)[totals$period]))
decomp <- readout(a38, "source_decomposition")
max_decomposition_error <- max(abs(decomp$total_difference - decomp$temporal_component - decomp$membership_component))
stopifnot(max_decomposition_error < 1e-12)
freq <- readout(a39, "inclusion_frequencies"); reps <- readout(a39, "simulation_replicates")
stopifnot(nrow(reps) == 24000, nrow(freq) == 2616, all(freq$inclusion_frequency >= 0 & freq$inclusion_frequency <= 1))
f <- freq %>% group_by(scenario, period) %>% summarise(n_selected = sum(inclusion_frequency), .groups = "drop")
stopifnot(all(abs(f$n_selected - 10) < 1e-10), all(reps$product_gap_pp >= -1e-10),
  all(reps$observed_list_gap_pp >= -1e-10), all(reps$day_gap_pp >= -1e-10), all(reps$level_gap_pp >= -1e-10))
# A citywide level offset and citywide activity multiplier cannot alter ranks.
base <- rowSums(E * A)
shift <- rowSums(10^((d$levels + 2) / 10) * (A * 1.3))
stopifnot(identical(order(-base), order(-shift)))
# Within-station persistent multiplicative errors preserve each station's hourly/product ratio.
scale <- exp(seq(-.4, .4, length.out = 109))
for (p in names(periods)) {
  h <- periods[[p]] + 1
  prod0 <- rowSums(A[, h, drop = FALSE]) * rowMeans(E[, h, drop = FALSE])
  joint0 <- rowSums((E * A)[, h, drop = FALSE])
  Ep <- E * scale; Ap <- A * rev(scale)
  prod1 <- rowSums(Ap[, h, drop = FALSE]) * rowMeans(Ep[, h, drop = FALSE])
  joint1 <- rowSums((Ep * Ap)[, h, drop = FALSE])
  stopifnot(max(abs(joint0 / prod0 - joint1 / prod1)) < 1e-10)
}
# Analyses 40-43: indicators, reference counts, profile cells, priority classes and date scenarios rebuilt independently.
a40 <- "analysis_40_assessment_benchmarks"; a41 <- "analysis_41_priority_profile_reference"
a42 <- "analysis_42_priority_places"; a43 <- "analysis_43_recording_dates"
ind <- readout(a40, "indicators"); ex <- readout(a40, "exceedance")
dbm <- function(h) 10 * log10(rowMeans(E[, h + 1, drop = FALSE]))
lden <- 10 * log10((12 * 10^(dbm(7:18) / 10) + 4 * 10^((dbm(19:22) + 5) / 10) + 8 * 10^((dbm(c(23, 0:6)) + 10) / 10)) / 24)
max_indicator_error <- max(abs(c(ind$lden - lden, ind$lday - dbm(7:22), ind$lnight - dbm(c(23, 0:6)),
  ind$lday_statutory - dbm(6:21), ind$lnight_statutory - dbm(c(22, 23, 0:5)))))
count_above <- function(ref, per, def) filter(ex, reference == ref, limit_period == per, definition == def)$n_above
stopifnot(max_indicator_error < 1e-10, count_above("GB 3096 class 4a", "night", "primary") == sum(dbm(c(23, 0:6)) > 55),
  count_above("GB 3096 class 4a", "night", "statutory") == sum(dbm(c(22, 23, 0:5)) > 55),
  count_above("GB 3096 class 4a", "day", "primary") == sum(dbm(7:22) > 70), count_above("WHO road traffic", "lden", "primary") == sum(lden > 53))
rules40 <- readout(a40, "rules")
stat_h <- list(day = 6:21 + 1, night = c(22, 23, 0:5) + 1); max_statutory_error <- 0
for (p in names(stat_h)) {
  target <- rowSums((E * A)[, stat_h[[p]], drop = FALSE]); ref_rank <- rank(-target, ties.method = "first")
  vals <- list(period_joint = target, full_joint = rowSums(E * A), period_level = rowMeans(E[, stat_h[[p]], drop = FALSE]),
    period_product = rowSums(A[, stat_h[[p]], drop = FALSE]) * rowMeans(E[, stat_h[[p]], drop = FALSE]), lden_level = lden)
  for (m in names(vals)) {
    z <- filter(rules40, definition == "statutory", period == p, method == m) %>% arrange(budget); r <- rank(-vals[[m]], ties.method = "first")
    stopifnot(nrow(z) == 30, all(z$overlap == sapply(1:30, function(k) sum(r <= k & ref_rank <= k))))
    max_statutory_error <- max(max_statutory_error, max(abs(z$fraction_index - sapply(1:30, function(k) sum(target[r <= k]) / sum(target)))))
  }
}
stopifnot(max_statutory_error < 1e-12)
al <- readout(a40, "alignment_ratio") %>% filter(definition == "raw50", period == "night") %>% arrange(station_key)
hn <- periods$night + 1
stopifnot(max(abs(al$ratio - rowSums((E * A)[, hn]) / (rowSums(A[, hn]) * rowMeans(E[, hn])))) < 1e-10)
prof <- readout(a41, "profile_by_size"); net <- readout(a41, "network_reference")
G <- d$groups; max_profile_error <- 0
for (k in c(1, 5, 10, 15, 20, 30)) {
  D <- rank(-rowSums((E * A)[, periods$day + 1]), ties.method = "first") <= k
  N <- rank(-rowSums((E * A)[, hn]), ties.method = "first") <= k
  z <- filter(prof, rule == "joint", budget == k, source == "traffic")
  cells <- c(mean(G[D, periods$day + 1] == "traffic"), mean(G[D, hn] == "traffic"), mean(G[N, periods$day + 1] == "traffic"), mean(G[N, hn] == "traffic"))
  max_profile_error <- max(max_profile_error, max(abs(cells - c(z$day_day, z$day_night, z$night_day, z$night_night))))
}
stopifnot(max_profile_error < 1e-12, max(abs(prof$total_difference - prof$temporal_component - prof$membership_component)) < 1e-12,
  abs(filter(net, period == "day", source == "traffic")$share - mean(G[, periods$day + 1] == "traffic")) < 1e-12,
  abs(filter(net, period == "night", source == "traffic")$share - mean(G[, hn] == "traffic")) < 1e-12)
rr <- readout(a41, "random_reference")
stopifnot(nrow(rr) == 120, all(rr$random_q025 <= rr$random_median), all(rr$random_median <= rr$random_q975))
pl <- readout(a42, "priority_places"); ng <- readout(a42, "night_gain_by_size")
g37 <- readout(a37, "gains") %>% filter(definition == "raw50", period == "night") %>% arrange(budget)
stopifnot(identical(as.integer(table(pl$priority_class)[c("day_only", "night_only", "other", "shared")]), c(4L, 4L, 95L, 6L)),
  sum(pl$night_louder) == sum(dbm(c(23, 0:6)) > dbm(7:22)), max(abs(ng$gain_over_full - g37$period_gain_joint)) < 1e-10,
  max(abs(ng$gain_over_day - g37$night_vs_day_list_gain)) < 1e-10)
od <- readout(a43, "omit_one_date"); de <- readout(a43, "date_error_replicates"); df <- readout(a43, "date_error_inclusion")
fs <- readout(a43, "full_sample"); oe <- readout(a43, "omit_one_entrance") %>% arrange(omitted)
o42 <- readout(a42, "night_gain_omit_one") %>% arrange(omitted)
stopifnot(identical(oe$omitted, o42$omitted), max(abs(oe$night_gain_over_full - o42$gain_over_full)) < 1e-10)
stopifnot(nrow(od) == 25, sum(od$omitted_entrances) == 109, all(od$entrances + od$omitted_entrances == 109),
  abs(fs$night_gain_over_full - g37$period_gain_joint[10]) < 1e-10, nrow(de) == 8000, all(de$product_gap_pp >= -1e-10),
  all(de$observed_list_gap_pp >= -1e-10), all(abs(tapply(df$inclusion_frequency, paste(df$scenario, df$period), sum) - 10) < 1e-10))
checks <- tibble(check = c("Independent score reconstruction", "Independent coverage reconstruction", "Source proportions and period lengths",
  "Exact source decomposition", "Simulation dimensions and non-negative objective gaps", "Ten selections in each simulation",
  "Common citywide error leaves ranks unchanged", "Persistent site errors preserve within-site alignment ratios",
  "Period indicators and reference-value counts", "Statutory-period rules reconstructed", "Source profile cells and network shares",
  "Priority classes and night gains agree with analysis 37", "Recording-date omissions and date-shared simulations"),
  status = "PASS", detail = c(paste("max relative error", signif(max_score_error, 3)),
    paste("max absolute error", signif(max_coverage_error, 3)), paste(nrow(totals), "site-period-label records"),
    paste("max error", signif(max_decomposition_error, 3)), "24,000 period simulations; 2,616 inclusion frequencies",
    "all 24 scenario-period totals equal 10", "analytical invariant verified", "all five periods verified",
    paste("max indicator error", signif(max_indicator_error, 3)), paste("max coverage error", signif(max_statutory_error, 3)),
    paste("max cell error", signif(max_profile_error, 3)), "6 shared, 4 day-only, 4 night-only; gains identical at sizes 1-30",
    "25 omissions; 8,000 period simulations; totals equal 10"))
write_csv(checks, file.path(root, a39, "outcome/integrated_verification.csv"))
print(checks, width = Inf)
