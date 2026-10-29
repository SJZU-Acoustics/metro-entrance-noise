suppressPackageStartupMessages(library(tidyverse))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- "intermediate"
source(file.path("code", "shared", "helpers.R"))
d <- load_p11(root)
out <- file.path(root, "analysis_40_assessment_benchmarks/outcome")
ind <- p11_indicators(d)
stopifnot(max(abs(ind$l24 - d$noise$laeq_24h_leq_db_a)) < 1e-8, all(ind$lden > ind$l24))
# Reference values on single-day outdoor records; not compliance determinations.
refs <- tribble(~reference, ~limit_period, ~value,
  "GB 3096 class 4a", "day", 70, "GB 3096 class 4a", "night", 55,
  "GB 3096 class 2", "day", 60, "GB 3096 class 2", "night", 50,
  "WHO road traffic", "night", 45, "WHO road traffic", "lden", 53)
exceedance <- pmap_dfr(refs, function(reference, limit_period, value) {
  cols <- if (limit_period == "lden") c(primary = "lden") else
    c(primary = paste0("l", limit_period), statutory = paste0("l", limit_period, "_statutory"))
  imap_dfr(cols, ~tibble(reference = reference, limit_period = limit_period, value_db = value, definition = .y,
    indicator = .x, n_entrances = nrow(ind), n_above = sum(ind[[.x]] > value),
    median_margin_db = median(ind[[.x]] - value), max_margin_db = max(ind[[.x]] - value)))
})
# Selection rules under each period definition, evaluated against that definition's joint index.
compare_rules <- function(periods, set, other_periods) {
  s <- p11_scores(d$panel, periods); s_other <- p11_scores(d$panel, other_periods)
  full <- filter(s, period == "full") %>% arrange(station_key)
  day <- filter(s, period == "day") %>% arrange(station_key)
  stopifnot(identical(full$station_key, ind$station_key))
  map_dfr(c("day", "night"), function(p) {
    target <- filter(s, period == p) %>% arrange(station_key)
    values <- list(full_level = full$level, period_level = target$level, lden_level = ind$lden,
      period_activity = target$raw50_activity, full_joint = full$raw50_joint, day_joint = day$raw50_joint,
      period_product = target$raw50_factorised, period_joint = target$raw50_joint,
      other_definition_joint = (filter(s_other, period == p) %>% arrange(station_key))$raw50_joint)
    ranks <- lapply(values, function(z) rank(-z, ties.method = "first"))
    index <- target$raw50_joint
    bind_rows(lapply(names(values), function(m) tibble(definition = set, period = p, method = m, budget = 1:30,
      overlap = map_int(1:30, ~sum(ranks[[m]] <= .x & ranks$period_joint <= .x)),
      fraction_index = map_dbl(1:30, ~sum(index[ranks[[m]] <= .x]) / sum(index))))) %>%
      mutate(rank_list = list(tibble(station_key = target$station_key, rank = ranks$period_joint)))
  })
}
primary_periods <- p11_periods()[c("full", "day", "night")]
rules_all <- bind_rows(compare_rules(primary_periods, "primary", p11_periods_statutory()),
  compare_rules(p11_periods_statutory(), "statutory", primary_periods))
rank_lists <- rules_all %>% distinct(definition, period, rank_list) %>% unnest(rank_list)
rules <- select(rules_all, -rank_list)
# The primary-period rules must reproduce analysis 37 exactly.
a37 <- read_csv(file.path(root, "analysis_37_information_choices/outcome/coverage.csv"), show_col_types = FALSE) %>%
  filter(definition == "raw50", period %in% c("day", "night")) %>% select(period, method, budget, a37_overlap = overlap, a37_fraction = fraction_index)
chk <- inner_join(filter(rules, definition == "primary"), a37, by = c("period", "method", "budget"))
stopifnot(nrow(chk) == 2 * 7 * 30, all(chk$overlap == chk$a37_overlap), max(abs(chk$fraction_index - chk$a37_fraction)) < 1e-12)
definition_agreement <- rank_lists %>% pivot_wider(names_from = definition, values_from = rank) %>%
  group_by(period) %>% summarise(budget = list(1:30), shared = list(map_int(1:30, ~sum(primary <= .x & statutory <= .x))),
    rho = cor(primary, statutory, method = "spearman"), .groups = "drop") %>% unnest(c(budget, shared))
# Hourly joint score over period product, by entrance: the only quantity through which alignment can reorder entrances.
s <- p11_scores(d$panel)
alignment <- map_dfr(c("raw50", "raw20", "sqrt50", "cap50"), function(v) s %>% filter(period != "full") %>%
  transmute(definition = v, period, station_key, joint = .data[[paste0(v, "_joint")]], product = .data[[paste0(v, "_factorised")]],
    ratio = joint / product))
alignment_summary <- alignment %>% group_by(definition, period) %>% summarise(ratio_median = median(ratio),
  ratio_min = min(ratio), ratio_max = max(ratio), sd_log_ratio = sd(log(ratio)), sd_log_product = sd(log(product)),
  rho_joint_product = cor(joint, product, method = "spearman"), .groups = "drop")
p11_save(list(indicators = ind, exceedance = exceedance, rules = rules, period_rank_lists = rank_lists,
  definition_agreement = definition_agreement, alignment_ratio = alignment, alignment_summary = alignment_summary), out, d, script, root)
print(exceedance, width = Inf)
print(rules %>% filter(budget == 10) %>% mutate(coverage = round(100 * fraction_index, 2)) %>% select(-fraction_index) %>%
  pivot_wider(names_from = c(definition, period), values_from = c(overlap, coverage)), width = Inf)
print(filter(definition_agreement, budget %in% c(5, 10, 20, 30)), width = Inf)
print(alignment_summary, n = Inf, width = Inf)
