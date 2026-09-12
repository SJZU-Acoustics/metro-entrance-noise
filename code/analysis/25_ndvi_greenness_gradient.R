# 25_ndvi_greenness_gradient.R
# Forward-looking analysis: satellite greenness (NDVI) and the soundscape.
#
# Question. The built-environment screening used land-use percentages, where
# parks/green land showed only weak associations. NDVI is a continuous
# remote-sensing greenness measure available at six buffer radii
# (20-500 m). Does surrounding greenness relate to lower noise and to a more
# natural source mix, and at which spatial scale is any greenness signal
# strongest? This speaks to green infrastructure as a quiet-soundscape lever.

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

# Spearman partial correlation of x and y given covariate z (rank residuals).
partial_spearman <- function(x, y, z) {
  ok <- complete.cases(x, y, z)
  if (sum(ok) < 8) return(NA_real_)
  rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
  ex <- residuals(lm(rx ~ rz))
  ey <- residuals(lm(ry ~ rz))
  suppressWarnings(cor(ex, ey))
}

master <- read_csv(file.path(data_dir, "station_master.csv"), show_col_types = FALSE)
noise <- read_csv(file.path(data_dir, "noise_metrics_station.csv"), show_col_types = FALSE) %>%
  select(station_key, laeq_24h_leq_db_a, lday_leq_db_a, lnight_leq_db_a, lpeak_leq_db_a)

dat <- master %>%
  select(-any_of("laeq_24h_leq_db_a")) %>%
  inner_join(noise, by = "station_key") %>%
  filter(!is.na(laeq_24h_leq_db_a))

ndvi_cols <- c("20_ndvi", "50_ndvi", "100_ndvi", "200_ndvi", "300_ndvi", "500_ndvi")
ndvi_scale <- c(20, 50, 100, 200, 300, 500)
outcomes <- c("laeq_24h_leq_db_a", "lday_leq_db_a", "lnight_leq_db_a", "lpeak_leq_db_a",
              "sound_natural_share", "sound_traffic_share")

# Spearman of each NDVI scale against each outcome.
rows <- list()
for (i in seq_along(ndvi_cols)) {
  nc <- ndvi_cols[i]
  for (oc in outcomes) {
    r <- cor_with_p(dat[[nc]], dat[[oc]], "spearman")
    rows[[paste(nc, oc)]] <- r %>%
      mutate(ndvi_scale_m = ndvi_scale[i], ndvi_col = nc, outcome = oc, .before = 1)
  }
}
ndvi_corr <- bind_rows(rows) %>%
  group_by(outcome) %>%
  mutate(q_fdr = p.adjust(p, method = "BH")) %>%
  ungroup()
write_csv(ndvi_corr, file.path(out_dir, "25_ndvi_gradient__correlations.csv"))

# Best (strongest |rho|) NDVI scale per outcome.
best_scale <- ndvi_corr %>%
  group_by(outcome) %>%
  slice_max(abs(rho), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(outcome)
write_csv(best_scale, file.path(out_dir, "25_ndvi_gradient__best_scale.csv"))

# Partial Spearman: NDVI vs noise/natural-share controlling for arterial road density.
ctrl <- "100m_urban_arterial_road_density_m_m"
prows <- list()
for (i in seq_along(ndvi_cols)) {
  nc <- ndvi_cols[i]
  for (oc in c("laeq_24h_leq_db_a", "lnight_leq_db_a", "sound_natural_share")) {
    prows[[paste(nc, oc)]] <- tibble(
      ndvi_scale_m = ndvi_scale[i], ndvi_col = nc, outcome = oc,
      partial_rho_ctrl_arterial = partial_spearman(dat[[nc]], dat[[oc]], dat[[ctrl]])
    )
  }
}
ndvi_partial <- bind_rows(prows)
write_csv(ndvi_partial, file.path(out_dir, "25_ndvi_gradient__partial_arterial.csv"))

# Bootstrap CI (1000) for the strongest NDVI-LAeq_24h association.
best_laeq <- ndvi_corr %>% filter(outcome == "laeq_24h_leq_db_a") %>%
  slice_max(abs(rho), n = 1, with_ties = FALSE)
bc <- best_laeq$ndvi_col[1]
set.seed(2026)
boot_n <- 1000
dd <- dat %>% filter(!is.na(.data[[bc]]), !is.na(laeq_24h_leq_db_a))
boot_rho <- replicate(boot_n, {
  idx <- sample.int(nrow(dd), replace = TRUE)
  safe_cor(dd[[bc]][idx], dd$laeq_24h_leq_db_a[idx], "spearman")
})
boot_ci <- tibble(
  association = paste0(bc, " ~ laeq_24h"),
  rho = safe_cor(dd[[bc]], dd$laeq_24h_leq_db_a, "spearman"),
  ci_lo = quantile(boot_rho, 0.025, na.rm = TRUE),
  ci_hi = quantile(boot_rho, 0.975, na.rm = TRUE),
  n_boot = boot_n
)
write_csv(boot_ci, file.path(out_dir, "25_ndvi_gradient__boot_ci.csv"))

# Figure a: rho across NDVI buffer scale for the main outcomes.
plot_outcomes <- c("laeq_24h_leq_db_a", "lnight_leq_db_a", "sound_natural_share")
oc_labels <- c(laeq_24h_leq_db_a = "LAeq,24h", lnight_leq_db_a = "Lnight",
               sound_natural_share = "Natural share")
grad_df <- ndvi_corr %>% filter(outcome %in% plot_outcomes) %>%
  mutate(outcome_lab = oc_labels[outcome])
p1 <- ggplot(grad_df, aes(ndvi_scale_m, rho, colour = outcome_lab, shape = outcome_lab)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("LAeq,24h" = "#0072B2", "Lnight" = "#D55E00",
                                 "Natural share" = "#009E73")) +
  scale_x_continuous(breaks = ndvi_scale) +
  labs(x = "NDVI buffer radius (m)", y = "Spearman rho") +
  theme_pub() +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))
ggsave(file.path(out_dir, "25_ndvi_gradient__rho_by_scale.png"),
       p1, width = 85 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

# Figure b: scatter of best-scale NDVI vs LAeq_24h.
p2 <- ggplot(dat, aes(.data[[bc]], laeq_24h_leq_db_a)) +
  geom_point(colour = "#0072B2", alpha = 0.7, size = 1.4) +
  geom_smooth(method = "loess", se = TRUE, colour = "#D55E00",
              fill = "#D55E00", alpha = 0.15, linewidth = 0.6) +
  labs(x = paste0("NDVI (", sub("_ndvi", " m", bc), ")"),
       y = expression(italic(L)[plain("Aeq,24h")]~"(dB)")) +
  theme_pub()
ggsave(file.path(out_dir, "25_ndvi_gradient__best_scatter.png"),
       p2, width = 85 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

message("Completed NDVI greenness gradient analysis; best LAeq scale = ", bc)
