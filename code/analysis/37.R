suppressPackageStartupMessages(library(tidyverse))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- "intermediate"
source(file.path("code", "shared", "helpers.R"))
d <- load_p11(root)
s <- p11_scores(d$panel)
out <- file.path(root, "analysis_37_information_choices/outcome")
methods <- c("full_level", "period_level", "full_joint", "period_joint", "period_product", "day_joint", "period_activity")
definitions <- c("raw50", "raw20", "sqrt50", "cap50")
targets <- c("day", "peak", "offpeak", "night")
rows <- list(); rankings <- list(); n <- 0L; j <- 0L
for (v in definitions) for (p in targets) {
  target <- filter(s, period == p) %>% arrange(station_key)
  full <- filter(s, period == "full") %>% arrange(station_key)
  day <- filter(s, period == "day") %>% arrange(station_key)
  stopifnot(identical(target$station_key, full$station_key), identical(target$station_key, day$station_key))
  values <- list(full_level = full$level, period_level = target$level,
    full_joint = full[[paste0(v, "_joint")]], period_joint = target[[paste0(v, "_joint")]],
    period_product = target[[paste0(v, "_factorised")]], day_joint = day[[paste0(v, "_joint")]],
    period_activity = target[[paste0(v, "_activity")]])
  ranks <- lapply(values, function(z) rank(-z, ties.method = "first"))
  target_index <- target[[paste0(v, "_joint")]]
  for (m in methods) {
    j <- j + 1L
    rankings[[j]] <- tibble(definition = v, period = p, method = m, station_key = target$station_key,
      station_name_cn = target$station_name_cn, rank = ranks[[m]], selection_score = values[[m]])
    for (k in 1:30) {
      n <- n + 1L; sel <- ranks[[m]] <= k; ref <- ranks$period_joint <= k
      rows[[n]] <- tibble(definition = v, period = p, method = m, budget = k,
        overlap = sum(sel & ref), fraction_index = sum(target_index[sel]) / sum(target_index),
        fraction_default_index = sum(target$raw50_joint[sel]) / sum(target$raw50_joint),
        fraction_activity = sum(target[[paste0(v, "_activity")]][sel]) / sum(target[[paste0(v, "_activity")]]),
        fraction_energy = sum(target$level[sel]) / sum(target$level))
    }
  }
}
coverage <- bind_rows(rows)
rankings <- bind_rows(rankings)
gains <- coverage %>% select(definition, period, budget, method, fraction_index) %>%
  pivot_wider(names_from = method, values_from = fraction_index) %>% mutate(
    period_gain_level = 100 * (period_level - full_level),
    period_gain_joint = 100 * (period_joint - full_joint),
    activity_gain_full = 100 * (full_joint - full_level),
    activity_gain_period = 100 * (period_joint - period_level),
    hourly_alignment_gain = 100 * (period_joint - period_product),
    night_vs_day_list_gain = 100 * (period_joint - day_joint))
agreement <- rankings %>% select(definition, period, method, station_key, rank) %>%
  pivot_wider(names_from = method, values_from = rank) %>% group_by(definition, period) %>% summarise(
    rho_hourly_product = cor(period_joint, period_product, method = "spearman"),
    overlap10_hourly_product = sum(period_joint <= 10 & period_product <= 10), .groups = "drop")
metadata <- tibble(item = c("n_sites", "n_hours", "cap95_unique_hours", "tie_rule"),
  value = c(nrow(d$master), nrow(d$panel), as.numeric(quantile(d$panel$pop_heat_50m, .95)), "station_key lexical order"))
p11_save(list(coverage = coverage, gains = gains, rankings = rankings, agreement = agreement, metadata = metadata), out, d, script, root)
print(gains %>% filter(definition == "raw50", budget == 10) %>%
  select(period, full_level, period_level, full_joint, period_joint, period_product, day_joint,
    period_gain_joint, activity_gain_period, hourly_alignment_gain), width = Inf)
print(agreement, n = Inf)
