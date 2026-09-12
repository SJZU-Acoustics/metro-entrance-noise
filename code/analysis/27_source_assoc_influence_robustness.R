# 27_source_assoc_influence_robustness.R
# Robustness analysis: are the headline source-composition associations stable?
#
# Question. The headline result is that traffic share is positively and natural
# share negatively associated with station noise. Both predictors are heavily
# skewed (traffic share median ~0.96; many zero-natural stations), so the
# associations could in principle be driven by a few stations or by the zero
# spike. This script stress-tests the headline correlations with alternative
# correlation methods, leave-one-out jackknife, sub-sample exclusions, and a
# bootstrap interval.

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
cor_p <- function(x, y, method) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 6) return(NA_real_)
  suppressWarnings(cor.test(x[ok], y[ok], method = method)$p.value)
}

master <- read_csv(file.path(data_dir, "station_master.csv"), show_col_types = FALSE)
noise <- read_csv(file.path(data_dir, "noise_metrics_station.csv"), show_col_types = FALSE) %>%
  select(station_key, laeq_24h_leq_db_a, lnight_leq_db_a)

dat <- master %>%
  select(station_key, station_name_cn, sound_traffic_share, sound_natural_share) %>%
  inner_join(noise, by = "station_key") %>%
  filter(!is.na(sound_traffic_share), !is.na(laeq_24h_leq_db_a))

# Define the headline associations to stress-test.
assoc <- tribble(
  ~label,                         ~xcol,                 ~ycol,
  "traffic_share ~ LAeq_24h",     "sound_traffic_share", "laeq_24h_leq_db_a",
  "natural_share ~ LAeq_24h",     "sound_natural_share", "laeq_24h_leq_db_a",
  "traffic_share ~ Lnight",       "sound_traffic_share", "lnight_leq_db_a",
  "natural_share ~ Lnight",       "sound_natural_share", "lnight_leq_db_a"
)

# 1. Method comparison: Pearson / Spearman / Kendall.
method_rows <- list()
for (i in seq_len(nrow(assoc))) {
  x <- dat[[assoc$xcol[i]]]; y <- dat[[assoc$ycol[i]]]
  for (m in c("pearson", "spearman", "kendall")) {
    method_rows[[paste(i, m)]] <- tibble(
      label = assoc$label[i], method = m,
      r = safe_cor(x, y, m), p = cor_p(x, y, m)
    )
  }
}
method_tab <- bind_rows(method_rows)
write_csv(method_tab, file.path(out_dir, "27_source_robustness__methods.csv"))

# 2. Leave-one-out jackknife on Spearman rho.
jack_summary <- list()
jack_long <- list()
for (i in seq_len(nrow(assoc))) {
  x <- dat[[assoc$xcol[i]]]; y <- dat[[assoc$ycol[i]]]
  ok <- complete.cases(x, y)
  xo <- x[ok]; yo <- y[ok]; keys <- dat$station_key[ok]
  full <- safe_cor(xo, yo, "spearman")
  loo <- vapply(seq_along(xo), function(j) safe_cor(xo[-j], yo[-j], "spearman"), numeric(1))
  loo_p <- vapply(seq_along(xo), function(j) cor_p(xo[-j], yo[-j], "spearman"), numeric(1))
  jack_long[[assoc$label[i]]] <- tibble(
    label = assoc$label[i], dropped_station = keys, loo_rho = loo, loo_p = loo_p
  )
  # Most influential = removal causing largest absolute change in rho.
  infl_idx <- which.max(abs(loo - full))
  jack_summary[[assoc$label[i]]] <- tibble(
    label = assoc$label[i],
    full_rho = full,
    loo_min = min(loo), loo_max = max(loo), loo_range = max(loo) - min(loo),
    n_sign_flip = sum(sign(loo) != sign(full)),
    n_lose_sig = sum(loo_p >= 0.05),
    most_influential_station = keys[infl_idx],
    most_influential_delta = loo[infl_idx] - full
  )
}
jack_summary_tab <- bind_rows(jack_summary)
write_csv(jack_summary_tab, file.path(out_dir, "27_source_robustness__jackknife_summary.csv"))
write_csv(bind_rows(jack_long), file.path(out_dir, "27_source_robustness__jackknife_loo.csv"))

# 3. Sub-sample exclusions.
sub_rows <- list()
for (i in seq_len(nrow(assoc))) {
  xc <- assoc$xcol[i]; yc <- assoc$ycol[i]
  full_d <- dat %>% filter(!is.na(.data[[xc]]), !is.na(.data[[yc]]))
  full <- safe_cor(full_d[[xc]], full_d[[yc]], "spearman")
  # Exclude zero-natural stations (gradient-only test).
  nz <- full_d %>% filter(sound_natural_share > 0)
  rho_nz <- safe_cor(nz[[xc]], nz[[yc]], "spearman")
  # Trim 5% most extreme traffic-share stations at each tail.
  q <- quantile(full_d$sound_traffic_share, c(0.05, 0.95), na.rm = TRUE)
  tr <- full_d %>% filter(sound_traffic_share >= q[1], sound_traffic_share <= q[2])
  rho_tr <- safe_cor(tr[[xc]], tr[[yc]], "spearman")
  sub_rows[[assoc$label[i]]] <- tibble(
    label = assoc$label[i],
    full_rho = full, n_full = nrow(full_d),
    rho_excl_zero_natural = rho_nz, n_excl_zero = nrow(nz),
    rho_trim5_traffic = rho_tr, n_trim = nrow(tr)
  )
}
sub_tab <- bind_rows(sub_rows)
write_csv(sub_tab, file.path(out_dir, "27_source_robustness__subsets.csv"))

# 4. Bootstrap CI (1000) for each headline Spearman rho.
set.seed(2026)
boot_n <- 1000
boot_rows <- list()
for (i in seq_len(nrow(assoc))) {
  d <- dat %>% filter(!is.na(.data[[assoc$xcol[i]]]), !is.na(.data[[assoc$ycol[i]]]))
  br <- replicate(boot_n, {
    idx <- sample.int(nrow(d), replace = TRUE)
    safe_cor(d[[assoc$xcol[i]]][idx], d[[assoc$ycol[i]]][idx], "spearman")
  })
  boot_rows[[assoc$label[i]]] <- tibble(
    label = assoc$label[i],
    rho = safe_cor(d[[assoc$xcol[i]]], d[[assoc$ycol[i]]], "spearman"),
    ci_lo = quantile(br, 0.025, na.rm = TRUE),
    ci_hi = quantile(br, 0.975, na.rm = TRUE),
    n_boot = boot_n
  )
}
boot_tab <- bind_rows(boot_rows)
write_csv(boot_tab, file.path(out_dir, "27_source_robustness__boot_ci.csv"))

# Figure: leave-one-out rho distributions for the two LAeq_24h headline pairs.
loo_plot_df <- bind_rows(jack_long) %>%
  filter(label %in% c("traffic_share ~ LAeq_24h", "natural_share ~ LAeq_24h"))
full_ref <- jack_summary_tab %>%
  filter(label %in% c("traffic_share ~ LAeq_24h", "natural_share ~ LAeq_24h")) %>%
  select(label, full_rho)
p1 <- ggplot(loo_plot_df, aes(loo_rho, fill = label)) +
  geom_histogram(binwidth = 0.005, colour = "white", linewidth = 0.1) +
  geom_vline(data = full_ref, aes(xintercept = full_rho, colour = label),
             linewidth = 0.5, linetype = "dashed", show.legend = FALSE) +
  scale_fill_manual(values = c("traffic_share ~ LAeq_24h" = "#D55E00",
                               "natural_share ~ LAeq_24h" = "#009E73")) +
  scale_colour_manual(values = c("traffic_share ~ LAeq_24h" = "#D55E00",
                                 "natural_share ~ LAeq_24h" = "#009E73")) +
  labs(x = "Leave-one-out Spearman rho", y = "Count") +
  theme_pub() +
  theme(legend.position = c(0.5, 0.98), legend.justification = c(0.5, 1))
ggsave(file.path(out_dir, "27_source_robustness__jackknife_distribution.png"),
       p1, width = 85 / 25.4, height = 85 / 25.4, dpi = 600, bg = "white")

message("Completed source-association robustness analysis.")
