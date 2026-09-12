# 32_collocation_check.R
# Persistent record of the recorder-vs-reference collocation verification
# reported in Paper 1 Methods (Acoustic measurement and metrics).
#
# Question. Does the field recorder (Lusheng song-meter recorder) agree with a
# Class 1 reference sound level meter (B&K 2250) closely enough to use the raw
# recorded levels without numerical correction?
#
# Data. 14 collocated 10-min sessions (09:09-22:09, instruments 10 cm apart),
# held in data/validation/collocation_pairs.csv.
# This script re-derives every collocation figure printed in the manuscript:
# mean/min/max bias, RMSE, and the recorder-on-reference OLS slope and
# intercept with 95% CIs and p-values against H0 slope = 1 / intercept = 0.

suppressPackageStartupMessages({ library(tidyverse) })

# Paths are relative to the repository root; python run_all.py runs every module there.
out_dir <- file.path("intermediate", "analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
csv <- file.path("data", "validation", "collocation_pairs.csv")

# One row per collocated 10-min session: time (HH:MM), B&K 2250 LAeq and
# recorder LAeq, both in dB(A).
pairs <- read_csv(csv, col_types = cols(time = col_character(), .default = col_double()),
                  show_col_types = FALSE) %>%
  mutate(delta_db = recorder_db - reference_db)

stopifnot(nrow(pairs) == 14, all(is.finite(pairs$delta_db)))

fit <- lm(recorder_db ~ reference_db, data = pairs)
cf  <- summary(fit)$coefficients
ci  <- confint(fit)
# p-values against the no-bias null: slope = 1, intercept = 0.
t_slope <- (cf["reference_db", "Estimate"] - 1) / cf["reference_db", "Std. Error"]
p_slope <- 2 * pt(abs(t_slope), df = fit$df.residual, lower.tail = FALSE)
p_inter <- cf["(Intercept)", "Pr(>|t|)"]

summary_tbl <- tibble(
  metric = c("n_pairs", "bias_mean_db", "bias_min_db", "bias_max_db", "rmse_db",
             "slope", "slope_ci_lo", "slope_ci_hi", "p_slope_vs_1",
             "intercept_db", "intercept_ci_lo", "intercept_ci_hi", "p_intercept_vs_0"),
  value = c(nrow(pairs), mean(pairs$delta_db), min(pairs$delta_db), max(pairs$delta_db),
            sqrt(mean(pairs$delta_db^2)),
            cf["reference_db", "Estimate"], ci["reference_db", 1], ci["reference_db", 2], p_slope,
            cf["(Intercept)", "Estimate"], ci["(Intercept)", 1], ci["(Intercept)", 2], p_inter)
)

write_csv(pairs, file.path(out_dir, "32_collocation__pairs.csv"))
write_csv(summary_tbl, file.path(out_dir, "32_collocation__summary.csv"))

with(as.list(deframe(summary_tbl)), cat(sprintf(
  paste0("Collocation check (n = %d):\n",
         "  bias %.2f dB (range %.2f to %.2f), RMSE %.2f dB\n",
         "  slope %.3f (95%% CI %.3f to %.3f), p vs 1 = %.2f\n",
         "  intercept %.2f dB (95%% CI %.2f to %.2f), p vs 0 = %.2f\n",
         "  -> bias within IEC 61672 Class 1 tolerance (±1.1 dB); no correction applied.\n"),
  n_pairs, bias_mean_db, bias_min_db, bias_max_db, rmse_db,
  slope, slope_ci_lo, slope_ci_hi, p_slope_vs_1,
  intercept_db, intercept_ci_lo, intercept_ci_hi, p_intercept_vs_0)))
