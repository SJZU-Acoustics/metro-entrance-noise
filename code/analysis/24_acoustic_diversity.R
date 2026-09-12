# 24_acoustic_diversity.R
# Forward-looking analysis: acoustic source diversity / soundscape descriptors.
#
# Question. Beyond loudness, how acoustically diverse is each metro-entrance
# soundscape, and is diversity associated with lower noise, greener and more
# mixed surroundings, and station function? Metro entrances are dominated by
# traffic, so source diversity is itself a soundscape-quality dimension that
# the level-only analyses do not capture.

suppressPackageStartupMessages({
  library(tidyverse)
})

# Paths are relative to the repository root; python run_all.py runs every module there.
data_dir <- file.path("intermediate", "clean")
out_dir <- file.path("intermediate", "analysis")
if (!dir.exists(data_dir)) stop("Run from the repository root after code/load_data.py has built intermediate/clean/", call. = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(0.10, "cm"),
      axis.title = element_text(size = axis_title_size, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      legend.text = element_text(size = base_size),
      legend.title = element_blank(),
      legend.key = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

safe_cor <- function(x, y, method) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 6 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

cor_with_p <- function(x, y, method = "spearman") {
  ok <- complete.cases(x, y)
  if (sum(ok) < 6 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  tibble(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# Diversity helpers on a vector of category shares (non-negative, summing ~1).
shannon_h <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}
gini_simpson <- function(p) {
  p <- p[is.finite(p)]
  if (sum(p) == 0) return(NA_real_)
  1 - sum(p^2)
}

master <- read_csv(file.path(data_dir, "station_master.csv"), show_col_types = FALSE)
noise <- read_csv(file.path(data_dir, "noise_metrics_station.csv"), show_col_types = FALSE) %>%
  select(station_key, laeq_24h_leq_db_a, lday_leq_db_a, lnight_leq_db_a, lpeak_leq_db_a)

share_cols <- c("sound_natural_share", "sound_human_share",
                "sound_traffic_share", "sound_mechanical_share")

div <- master %>%
  select(-laeq_24h_leq_db_a) %>%
  inner_join(noise, by = "station_key") %>%
  filter(!is.na(laeq_24h_leq_db_a)) %>%
  rowwise() %>%
  mutate(
    shannon_diversity = shannon_h(c_across(all_of(share_cols))),
    gini_simpson = gini_simpson(c_across(all_of(share_cols))),
    source_richness = sum(c_across(all_of(share_cols)) > 0, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    pielou_evenness = ifelse(source_richness > 1,
                             shannon_diversity / log(source_richness), NA_real_),
    event_richness = active_event_type_count,
    dominance = dominant_event_share,
    inv_dominance = 1 - dominant_event_share
  )

div_metrics <- div %>%
  select(station_key, station_name_cn, line_name, station_type_code,
         station_type_name_en, laeq_24h_leq_db_a, lday_leq_db_a, lnight_leq_db_a,
         sound_traffic_share, sound_natural_share,
         shannon_diversity, gini_simpson, source_richness, pielou_evenness,
         event_richness, dominance, inv_dominance,
         ndvi_500 = `500_ndvi`, parks_green_500 = `500m_parks_green_space_land_pct`,
         shannon_landuse_300 = `300m_shannon_entropy_index`,
         arterial_100 = `100m_urban_arterial_road_density_m_m`)
write_csv(div_metrics, file.path(out_dir, "24_acoustic_diversity__station_metrics.csv"))

# Distribution summary of the diversity descriptors.
div_summary <- div %>%
  summarise(
    n_stations = n(),
    shannon_mean = mean(shannon_diversity, na.rm = TRUE),
    shannon_median = median(shannon_diversity, na.rm = TRUE),
    shannon_p25 = quantile(shannon_diversity, 0.25, na.rm = TRUE),
    shannon_p75 = quantile(shannon_diversity, 0.75, na.rm = TRUE),
    gini_median = median(gini_simpson, na.rm = TRUE),
    pielou_median = median(pielou_evenness, na.rm = TRUE),
    event_richness_median = median(event_richness, na.rm = TRUE),
    dominance_median = median(dominance, na.rm = TRUE),
    n_monotone = sum(shannon_diversity < 0.10, na.rm = TRUE),
    share_monotone = mean(shannon_diversity < 0.10, na.rm = TRUE),
    n_single_source = sum(source_richness == 1, na.rm = TRUE)
  )
write_csv(div_summary, file.path(out_dir, "24_acoustic_diversity__summary.csv"))

# Diversity vs acoustic level and context, Spearman with p.
div_vars <- c("shannon_diversity", "gini_simpson", "pielou_evenness",
              "event_richness", "inv_dominance")
ctx_vars <- c("laeq_24h_leq_db_a", "lday_leq_db_a", "lnight_leq_db_a",
              "ndvi_500", "parks_green_500", "shannon_landuse_300", "arterial_100")

corr_rows <- list()
for (dv in div_vars) {
  for (cv in ctx_vars) {
    r <- cor_with_p(div_metrics[[dv]], div_metrics[[cv]], "spearman")
    corr_rows[[paste(dv, cv)]] <- r %>% mutate(diversity_metric = dv, context = cv, .before = 1)
  }
}
corr_tab <- bind_rows(corr_rows) %>%
  group_by(diversity_metric) %>%
  mutate(q_fdr = p.adjust(p, method = "BH")) %>%
  ungroup()
write_csv(corr_tab, file.path(out_dir, "24_acoustic_diversity__correlations.csv"))

# Group differences in Shannon diversity by station type and by line.
kw_type <- with(div_metrics, kruskal.test(shannon_diversity ~ factor(station_type_code)))
kw_line <- with(div_metrics, kruskal.test(shannon_diversity ~ factor(line_name)))
group_tab <- tibble(
  grouping = c("station_type", "line"),
  statistic = c(unname(kw_type$statistic), unname(kw_line$statistic)),
  df = c(unname(kw_type$parameter), unname(kw_line$parameter)),
  p = c(kw_type$p.value, kw_line$p.value)
)
write_csv(group_tab, file.path(out_dir, "24_acoustic_diversity__group_tests.csv"))

# Bootstrap CI (1000) for the headline Shannon-LAeq_24h Spearman correlation.
set.seed(2026)
boot_n <- 1000
dd <- div_metrics %>% filter(!is.na(shannon_diversity), !is.na(laeq_24h_leq_db_a))
boot_rho <- replicate(boot_n, {
  idx <- sample.int(nrow(dd), replace = TRUE)
  safe_cor(dd$shannon_diversity[idx], dd$laeq_24h_leq_db_a[idx], "spearman")
})
boot_ci <- tibble(
  association = "shannon_diversity ~ laeq_24h",
  rho = safe_cor(dd$shannon_diversity, dd$laeq_24h_leq_db_a, "spearman"),
  ci_lo = quantile(boot_rho, 0.025, na.rm = TRUE),
  ci_hi = quantile(boot_rho, 0.975, na.rm = TRUE),
  n_boot = boot_n
)
write_csv(boot_ci, file.path(out_dir, "24_acoustic_diversity__boot_ci.csv"))

# Figure a: Shannon diversity vs LAeq_24h.
p1 <- ggplot(div_metrics, aes(laeq_24h_leq_db_a, shannon_diversity)) +
  geom_point(colour = "#0072B2", alpha = 0.7, size = 1.4) +
  geom_smooth(method = "loess", se = TRUE, colour = "#D55E00",
              fill = "#D55E00", alpha = 0.15, linewidth = 0.6) +
  labs(x = expression(italic(L)[plain("Aeq,24h")]~"(dB)"), y = "Source Shannon diversity (nats)") +
  theme_pub()
ggsave(file.path(out_dir, "24_acoustic_diversity__shannon_vs_laeq.png"),
       p1, width = 85 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

# Figure b: distribution of source richness.
rich_df <- div_metrics %>% count(source_richness, name = "n_stations")
p2 <- ggplot(rich_df, aes(factor(source_richness), n_stations)) +
  geom_col(fill = "#009E73", width = 0.7) +
  labs(x = "Distinct source categories (n)", y = "Stations (n)") +
  theme_pub()
ggsave(file.path(out_dir, "24_acoustic_diversity__source_richness.png"),
       p2, width = 85 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

message("Completed acoustic diversity analysis for ", nrow(div_metrics), " stations.")
