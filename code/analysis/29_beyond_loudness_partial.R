# 29_beyond_loudness_partial.R
# Decision test for whether a second (soundscape-quality) paper is independent
# of the first (level/burden) paper.
#
# Question. Do greenness (NDVI 500 m) and road density (arterial 100 m) still
# explain soundscape QUALITY -- source diversity and natural-source share --
# after controlling for LOUDNESS (LAeq_24h, Lnight)? If the quality signal
# survives level adjustment, Paper 2 has a thesis Paper 1 cannot make. If it
# collapses, "quality" is just inverse loudness and belongs inside Paper 1.

suppressPackageStartupMessages({ library(tidyverse) })

# Paths are relative to the repository root; python run_all.py runs every module there.
out_dir <- file.path("intermediate", "analysis")
if (!dir.exists(out_dir)) stop("Run from the repository root after the analysis modules", call. = FALSE)

dat <- read_csv(file.path(out_dir, "24_acoustic_diversity__station_metrics.csv"),
                show_col_types = FALSE)

# Partial Spearman via rank residualisation (generalises to >1 control).
# rank-then-Pearson on residuals == partial Spearman; p from t with df = n-2-k.
partial_spearman <- function(d, x, y, z = character(0)) {
  vars <- c(x, y, z)
  dd <- d[, vars, drop = FALSE]
  dd <- dd[complete.cases(dd), , drop = FALSE]
  n <- nrow(dd); k <- length(z)
  dr <- as.data.frame(lapply(dd, rank)); names(dr) <- vars
  if (k == 0) {
    r <- cor(dr[[x]], dr[[y]])
  } else {
    rx <- resid(lm(reformulate(z, response = x), data = dr))
    ry <- resid(lm(reformulate(z, response = y), data = dr))
    r <- cor(rx, ry)
  }
  dfree <- n - 2 - k
  tval <- r * sqrt(dfree / (1 - r^2))
  tibble(rho = r, p = 2 * pt(-abs(tval), dfree), n = n, k = k)
}

run_block <- function(d, label) {
  outcomes  <- c("sound_natural_share", "shannon_diversity", "gini_simpson", "inv_dominance")
  predictors <- c("ndvi_500", "arterial_100")
  rows <- list()
  for (y in outcomes) for (x in predictors) {
    # control sets: none, level only (two ways), level + the OTHER context driver
    other <- setdiff(predictors, x)
    ctrls <- list(
      none              = character(0),
      laeq              = "laeq_24h_leq_db_a",
      lnight            = "lnight_leq_db_a",
      laeq_plus_other   = c("laeq_24h_leq_db_a", other),
      lnight_plus_other = c("lnight_leq_db_a", other)
    )
    for (cn in names(ctrls)) {
      r <- partial_spearman(d, x, y, ctrls[[cn]])
      rows[[paste(y, x, cn)]] <- r %>%
        mutate(subset = label, outcome = y, predictor = x, control = cn, .before = 1)
    }
  }
  bind_rows(rows)
}

full <- run_block(dat, "all_109")
# Sensitivity: drop the zero-inflation in natural sound (50/109 monotone).
nz   <- run_block(dat %>% filter(sound_natural_share > 0), "natural_share_gt0")

res <- bind_rows(full, nz) %>%
  mutate(across(c(rho, p), ~round(., 3)))
write_csv(res, file.path(out_dir, "29_beyond_loudness_partial.csv"))

# Bootstrap CI (1000) for the two headline partials, all 109 stations.
set.seed(2026)
boot_partial <- function(d, x, y, z, B = 1000) {
  est <- partial_spearman(d, x, y, z)$rho
  bs <- replicate(B, {
    idx <- sample.int(nrow(d), replace = TRUE)
    partial_spearman(d[idx, ], x, y, z)$rho
  })
  tibble(association = paste0(y, " ~ ", x, " | ", paste(z, collapse = "+")),
         rho = round(est, 3),
         ci_lo = round(quantile(bs, .025, na.rm = TRUE), 3),
         ci_hi = round(quantile(bs, .975, na.rm = TRUE), 3), n_boot = B)
}
boot <- bind_rows(
  boot_partial(dat, "ndvi_500", "sound_natural_share", "laeq_24h_leq_db_a"),
  boot_partial(dat, "ndvi_500", "shannon_diversity",   "laeq_24h_leq_db_a"),
  boot_partial(dat, "ndvi_500", "sound_natural_share", c("laeq_24h_leq_db_a", "arterial_100")),
  boot_partial(dat, "ndvi_500", "shannon_diversity",   c("laeq_24h_leq_db_a", "arterial_100"))
)
write_csv(boot, file.path(out_dir, "29_beyond_loudness_partial__boot_ci.csv"))

# ---- console verdict ----------------------------------------------------
cat("\n==== ZERO-ORDER vs PARTIAL (all 109 stations) ====\n")
full %>%
  filter(control %in% c("none", "laeq", "lnight", "laeq_plus_other")) %>%
  transmute(outcome, predictor, control,
            rho = round(rho, 3), p = round(p, 3), n) %>%
  pivot_wider(names_from = control, values_from = c(rho, p)) %>%
  select(outcome, predictor,
         rho_none, rho_laeq, rho_lnight, rho_laeq_plus_other,
         p_none, p_laeq, p_lnight, p_laeq_plus_other) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n==== SENSITIVITY: natural_share > 0 only ====\n")
nz %>%
  filter(outcome %in% c("sound_natural_share", "shannon_diversity"),
         control %in% c("none", "laeq")) %>%
  transmute(outcome, predictor, control, rho = round(rho, 3),
            p = round(p, 3), n) %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n==== BOOTSTRAP CI for headline partials (all 109) ====\n")
as.data.frame(boot) %>% print(row.names = FALSE)
cat("\nWrote 29_beyond_loudness_partial.csv and __boot_ci.csv\n")
