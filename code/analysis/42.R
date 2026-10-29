suppressPackageStartupMessages(library(tidyverse))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- "intermediate"
source(file.path("code", "shared", "helpers.R"))
d <- load_p11(root)
out <- file.path(root, "analysis_42_priority_places/outcome")
ind <- p11_indicators(d); n <- nrow(d$master)
s <- p11_scores(d$panel) %>% arrange(period, station_key)
get <- function(p, col) filter(s, period == p)[[col]]
rk <- function(v) rank(-v, ties.method = "first")
r_day <- rk(get("day", "raw50_joint")); r_night <- rk(get("night", "raw50_joint")); r_full <- rk(get("full", "raw50_joint"))
day_h <- p11_periods()$day + 1; night_h <- p11_periods()$night + 1
grp <- function(h, g) rowMeans(d$groups[, h, drop = FALSE] == g)
places <- d$master %>% transmute(station_key, station_name_cn, line_name, station_type_name_en, measurement_date,
  day_type, lng_round, lat_round) %>%
  mutate(rank_day = r_day, rank_night = r_night, rank_full = r_full,
    priority_class = case_when(rank_day <= 10 & rank_night <= 10 ~ "shared", rank_day <= 10 ~ "day_only",
      rank_night <= 10 ~ "night_only", TRUE ~ "other"),
    lday = ind$lday, lnight = ind$lnight, night_minus_day = lnight - lday, night_louder = night_minus_day > 0,
    activity_day = get("day", "raw50_activity"), activity_night = get("night", "raw50_activity"),
    night_index_share = get("night", "raw50_joint") / sum(get("night", "raw50_joint")),
    day_index_share = get("day", "raw50_joint") / sum(get("day", "raw50_joint")),
    traffic_day = grp(day_h, "traffic"), human_day = grp(day_h, "human"), traffic_night = grp(night_h, "traffic"),
    natural_night = grp(night_h, "natural")) %>% bind_cols(as_tibble(d$env))
stopifnot(identical(places$station_key, ind$station_key), all(table(places$priority_class)[c("shared", "day_only", "night_only")] == c(6, 4, 4)))
vars <- c("lday", "lnight", "night_minus_day", "activity_day", "activity_night", "human_day", "traffic_night", d$env_names)
class_summary <- places %>% select(priority_class, all_of(vars)) %>% pivot_longer(-priority_class, names_to = "variable") %>%
  group_by(priority_class, variable) %>% summarise(entrances = n(), median = median(value), minimum = min(value), maximum = max(value), .groups = "drop")
# Night-specific list against lists reused from the 24-hour or daytime index, at every size.
idx <- get("night", "raw50_joint"); covf <- function(r, k, index = idx) 100 * sum(index[r <= k]) / sum(index)
night_gain <- tibble(budget = 1:30) %>% mutate(night_coverage = map_dbl(budget, ~covf(r_night, .x)),
  full_list_coverage = map_dbl(budget, ~covf(r_full, .x)), day_list_coverage = map_dbl(budget, ~covf(r_day, .x)),
  gain_over_full = night_coverage - full_list_coverage, gain_over_day = night_coverage - day_list_coverage,
  shared_with_full = map_int(budget, ~sum(r_night <= .x & r_full <= .x)), shared_with_day = map_int(budget, ~sum(r_night <= .x & r_day <= .x)))
# Omit each entrance belonging to any of the three ten-entrance lists, then reselect among the remaining 108.
listed <- which(r_day <= 10 | r_night <= 10 | r_full <= 10)
omit_one <- map_dfr(listed, function(i) {
  keep <- setdiff(seq_len(n), i); rn <- rk(idx[keep]); rf <- rk(get("full", "raw50_joint")[keep]); rd <- rk(get("day", "raw50_joint")[keep])
  cv <- function(r) 100 * sum(idx[keep][r <= 10]) / sum(idx[keep])
  tibble(omitted = d$master$station_key[i], gain_over_full = cv(rn) - cv(rf), gain_over_day = cv(rn) - cv(rd),
    shared_with_full = sum(rn <= 10 & rf <= 10), shared_with_day = sum(rn <= 10 & rd <= 10))
})
p11_save(list(priority_places = places, class_summary = class_summary, night_gain_by_size = night_gain,
  night_gain_omit_one = omit_one), out, d, script, root)
print(places %>% filter(priority_class != "other" | night_louder) %>% arrange(priority_class, rank_day) %>%
  transmute(station_key, station_name_cn, measurement_date, priority_class, rank_day, rank_night, lday = round(lday, 1), lnight = round(lnight, 1),
    diff = round(night_minus_day, 1), night_share = round(100 * night_index_share, 2)), n = Inf, width = Inf)
cat("night louder:", sum(places$night_louder), "of", n, "; median night minus day", round(median(places$night_minus_day), 2), "\n")
print(class_summary %>% mutate(across(median:maximum, ~signif(.x, 4))) %>% select(-entrances) %>%
  pivot_wider(names_from = priority_class, values_from = c(median, minimum, maximum)) %>% select(variable, starts_with("median")), n = Inf, width = Inf)
print(filter(night_gain, budget %in% c(5, 8, 9, 10, 11, 12, 15, 20, 30)) %>% mutate(across(where(is.double), ~round(.x, 2))), width = Inf)
print(summarise(omit_one, across(gain_over_full:shared_with_day, list(min = min, max = max))), width = Inf)
