# Shared analytical definitions and workbook-derived data entry.
suppressPackageStartupMessages(library(tidyverse))

load_p11 <- function(exploration_root) {
  paths <- file.path(exploration_root, "data/clean", c(
    "station_master.csv", "station_hour_panel.csv", "noise_metrics_station.csv"))
  master <- read_csv(paths[1], show_col_types = FALSE)
  panel <- read_csv(paths[2], show_col_types = FALSE)
  noise <- read_csv(paths[3], show_col_types = FALSE)
  good <- panel %>% group_by(station_key) %>%
    summarise(ok = n() == 24 && n_distinct(hour) == 24 &&
      all(is.finite(laeq_hour_db_a)) && all(!is.na(primary_event_cn)), .groups = "drop") %>%
    filter(ok) %>% pull(station_key)
  master <- master %>% filter(station_key %in% good) %>% arrange(station_key)
  panel <- panel %>% filter(station_key %in% good) %>% arrange(station_key, hour)
  noise <- noise %>% filter(station_key %in% good) %>% arrange(station_key)
  source_map <- c("鸟鸣声" = "natural", "虫鸣声" = "natural",
    "脚步声" = "human", "交谈声" = "human", "音乐声" = "human",
    "叫卖声" = "human", "店铺喇叭声" = "human", "汽车行驶声" = "traffic",
    "汽车鸣笛声" = "traffic", "摩托车声" = "traffic", "刹车声" = "traffic",
    "打扫声" = "mechanical", "机械摩擦声" = "mechanical")
  panel$source_group <- unname(source_map[panel$primary_event_cn])
  panel$secondary_group <- unname(source_map[panel$secondary_event_cn])
  panel$early_category_version <- panel$station_key %in% paste0("M_", c(93:99, 101))
  stopifnot(nrow(master) == 109, nrow(panel) == 109 * 24,
    !anyNA(panel$source_group), identical(master$station_key, noise$station_key))
  master$poi_total_density_300m <- rowSums(master[, endsWith(names(master), "_poi_density_300m")])
  env_names <- c("100m_urban_arterial_road_density_m_m", "500m_parks_green_space_land_pct",
    "500m_residential_land_pct", "500m_commercial_services_land_pct",
    "300m_parks_green_space_land_pct", "100m_shannon_entropy_index",
    "poi_total_density_300m", "500_ndvi")
  env <- as.matrix(master[, env_names])
  stopifnot(all(is.finite(env)))
  sources <- c("traffic", "natural", "human", "mechanical")
  source_matrix <- sapply(sources, function(g) rowMeans(matrix(panel$source_group == g,
    nrow = nrow(master), byrow = TRUE)))
  old_source <- as.matrix(master[, paste0("sound_", sources, "_share")])
  stopifnot(max(abs(source_matrix - old_source)) < 1e-10)
  level <- matrix(panel$laeq_hour_db_a, nrow = nrow(master), byrow = TRUE)
  group <- matrix(panel$source_group, nrow = nrow(master), byrow = TRUE)
  leq <- function(x) 10 * log10(rowMeans(10^(x / 10)))
  stopifnot(max(abs(leq(level) - noise$laeq_24h_leq_db_a)) < 1e-8)
  hashes <- tibble(file = basename(paths), sha256 = vapply(paths,
    function(p) digest::digest(file = p, algo = "sha256"), character(1)))
  list(master = master, panel = panel, noise = noise, env = env, env_names = env_names,
    source = source_matrix, groups = group, levels = level, sources = sources,
    input_hashes = hashes, leq = leq)
}

p11_periods <- function() list(full = 0:23, day = 7:22,
  peak = c(7, 8, 17, 18), offpeak = 10:16, night = c(23, 0:6))

# Default statutory day and night in China (06:00-22:00, 22:00-06:00), used as a sensitivity definition.
p11_periods_statutory <- function() list(full = 0:23, day = 6:21, night = c(22, 23, 0:5))

# Station-level equivalent levels under both period definitions, and the European day-evening-night level.
p11_indicators <- function(d) {
  e <- 10^(d$levels / 10); db <- function(h) 10 * log10(rowMeans(e[, h + 1, drop = FALSE]))
  eu <- p11_periods(); cn <- p11_periods_statutory()
  lden <- 10 * log10((12 * rowMeans(e[, 7:18 + 1]) + 4 * 10^.5 * rowMeans(e[, 19:22 + 1]) +
    8 * 10 * rowMeans(e[, c(23, 0:6) + 1])) / 24)
  tibble(station_key = d$master$station_key, l24 = db(0:23), lday = db(eu$day), lnight = db(eu$night),
    lday_statutory = db(cn$day), lnight_statutory = db(cn$night), lden = lden)
}

# Period definitions and activity transforms for the integrated campaign.
# The cap is fixed once on the 2,616 unique station-hours, not on duplicated periods.
p11_scores <- function(panel, periods = p11_periods()) {
  stopifnot(all(is.finite(panel$pop_heat_50m)), all(panel$pop_heat_50m >= 0),
    all(is.finite(panel$pop_heat_20m)), all(panel$pop_heat_20m >= 0))
  cap <- as.numeric(quantile(panel$pop_heat_50m, .95))
  x <- panel %>% mutate(energy = 10^(laeq_hour_db_a / 10),
    raw50 = pop_heat_50m, raw20 = pop_heat_20m, sqrt50 = sqrt(pop_heat_50m),
    cap50 = pmin(pop_heat_50m, cap))
  x <- imap_dfr(periods, ~filter(x, hour %in% .x) %>% mutate(period = .y))
  x %>% group_by(station_key, station_name_cn, period) %>% summarise(n_hours = n(),
    level = mean(energy), across(c(raw50, raw20, sqrt50, cap50),
      list(activity = ~sum(.x), joint = ~sum(.x * energy),
        factorised = ~sum(.x) * mean(energy)), .names = "{.col}_{.fn}"),
    .groups = "drop") %>% arrange(period, station_key)
}

p11_save <- function(outputs, out, d, script, root) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  iwalk(outputs, ~write_csv(.x, file.path(out, paste0(.y, ".csv"))))
  extra <- c(normalizePath(script), file.path("code", "shared", "helpers.R"))
  hashes <- bind_rows(d$input_hashes, tibble(file = sub(paste0(root, "/"), "", extra, fixed = TRUE),
    sha256 = vapply(extra, function(p) digest::digest(file = p, algo = "sha256"), character(1))))
  write_csv(hashes, file.path(out, "input_sha256.csv"))
  writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
}
