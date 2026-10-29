suppressPackageStartupMessages(library(tidyverse))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- "intermediate"
source(file.path("code", "shared", "helpers.R"))
d <- load_p11(root)
out <- file.path(root, "analysis_38_priority_sources/outcome")
profiles <- list(); decompositions <- list(); inventories <- list(); details <- list(); concentrations <- list()
for (sample in c("all109", "later101")) {
  panel <- if (sample == "all109") d$panel else filter(d$panel, !early_category_version)
  s <- p11_scores(panel) %>% group_by(period) %>% mutate(priority_rank = rank(-raw50_joint, ties.method = "first"),
    level_rank = rank(-level, ties.method = "first")) %>% ungroup()
  period_panel <- imap_dfr(p11_periods(), ~filter(panel, hour %in% .x) %>% mutate(period = .y))
  group_long <- period_panel %>% select(station_key, station_name_cn, period, hour, source_group, secondary_group) %>%
    pivot_longer(c(source_group, secondary_group), names_to = "label", values_to = "source")
  stopifnot(!anyNA(group_long$source))
  station_profiles <- group_long %>% count(station_key, station_name_cn, period, label, source, name = "hours") %>%
    group_by(station_key, station_name_cn, period, label) %>% complete(source = d$sources, fill = list(hours = 0)) %>%
    mutate(share = hours / sum(hours)) %>% ungroup() %>% left_join(s %>% select(station_key, period, priority_rank, level_rank), by = c("station_key", "period"))
  profiles[[sample]] <- station_profiles %>% mutate(sample = sample)
  for (k in c(5, 10, 20)) {
    day_ids <- filter(s, period == "day", priority_rank <= k)$station_key
    night_ids <- filter(s, period == "night", priority_rank <= k)$station_key
    # Four cells separate selected membership from the time at which the same sites are described.
    cells <- bind_rows(filter(station_profiles, station_key %in% day_ids, period %in% c("day", "night")) %>% mutate(selected_set = "day"),
      filter(station_profiles, station_key %in% night_ids, period %in% c("day", "night")) %>% mutate(selected_set = "night")) %>%
      group_by(label, source, selected_set, period) %>% summarise(share = mean(share), .groups = "drop")
    wide <- cells %>% mutate(cell = paste(selected_set, period, sep = "_")) %>% select(-selected_set, -period) %>%
      pivot_wider(names_from = cell, values_from = share) %>% mutate(sample = sample, budget = k,
        total_difference = night_night - day_day,
        temporal_component = .5 * ((day_night - day_day) + (night_night - night_day)),
        membership_component = .5 * ((night_day - day_day) + (night_night - day_night)))
    stopifnot(max(abs(wide$total_difference - wide$temporal_component - wide$membership_component)) < 1e-12)
    decompositions[[paste(sample, k)]] <- wide
  }
  selected <- filter(station_profiles, priority_rank <= 10, period != "full")
  agg <- selected %>% group_by(period, label, source) %>% summarise(n_sites = n(), total_hours = sum(hours),
    mean_share = mean(share), min_share = min(share), max_share = max(share),
    sites_present = sum(hours > 0), largest_site_fraction = if_else(sum(hours) > 0, max(hours) / sum(hours), NA_real_),
    loo_min = min((sum(share) - share) / (n() - 1)), loo_max = max((sum(share) - share) / (n() - 1)),
    .groups = "drop") %>% mutate(sample = sample)
  concentrations[[sample]] <- agg
  ids <- s %>% filter(period %in% c("day", "night"), priority_rank <= 10 | level_rank <= 10) %>%
    group_by(station_key) %>% summarise(selection_reason = paste(paste(period, if_else(priority_rank <= 10, "joint", ""),
      if_else(level_rank <= 10, "level", "")), collapse = "; "), .groups = "drop")
  inventories[[sample]] <- panel %>% inner_join(ids, by = "station_key") %>% mutate(sample = sample) %>%
    select(sample, station_key, station_name_cn, measurement_date, hour, primary_event_cn, secondary_event_cn,
      source_group, secondary_group, early_category_version, selection_reason)
  details[[sample]] <- period_panel %>% inner_join(filter(s, period != "full", priority_rank <= 10) %>%
    select(station_key, period, priority_rank), by = c("station_key", "period")) %>% mutate(sample = sample) %>%
    select(sample, period, station_key, station_name_cn, hour, primary_event_cn, secondary_event_cn, priority_rank)
}
profiles <- bind_rows(profiles); decomposition <- bind_rows(decompositions); concentration <- bind_rows(concentrations)
priority_profile <- profiles %>% filter(priority_rank <= 10, label == "source_group", period != "full") %>%
  group_by(sample, period, station_key, station_name_cn, priority_rank) %>%
  summarise(largest_group_share = max(share), tied_largest = sum(share == max(share)),
    largest_group = paste(source[share == max(share)], collapse = ";"), minority_hours = sum(hours) - max(hours), .groups = "drop")
p11_save(list(station_profiles = profiles, source_decomposition = decomposition,
  source_summary = concentration, priority_profile = priority_profile,
  priority_hours = bind_rows(details), verification_inventory = bind_rows(inventories)), out, d, script, root)
print(concentration %>% filter(sample == "all109", label == "source_group", period %in% c("day", "night")), width = Inf)
print(decomposition %>% filter(sample == "all109", budget == 10, label == "source_group"), width = Inf)
