# =============================================================================
# si_tables.R -- Generate the SI table fragments from the evidence bank.
# Writes BARE table bodies (tabular / longtable rows, no float, no caption,
# no label) to tables/si/, named by analysis; the SI .tex owns every caption.
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse) })
# Paths are relative to the repository root; run_all.py runs this script there.
EXPLORE_OUT   <- file.path("intermediate", "analysis")
EXPLORE_CLEAN <- file.path("intermediate", "clean")
TBL_DIR       <- file.path("output", "tables", "si")
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3)
sgn <- function(x) ifelse(x >= 0, paste0("$+$", f2(x)), paste0("$-$", f2(abs(x))))
pfmt <- function(p) ifelse(p < 0.001, "$<$0.001", formatC(p, format = "f", digits = 3))
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
wr <- function(lines, file) { writeLines(lines, file.path(TBL_DIR, file));
  message("  wrote ", file) }

# -----------------------------------------------------------------------------
# Station-level acoustic-metric summary (12 metrics)
# -----------------------------------------------------------------------------
nm <- read_csv(file.path(EXPLORE_CLEAN, "noise_metrics_station.csv"), show_col_types = FALSE) %>%
  filter(!is.na(laeq_24h_leq_db_a))
metrics <- tribble(
  ~col,                         ~lab,                  ~unit,
  "laeq_24h_leq_db_a",          "$L_\\mathrm{Aeq,24h}$","dB(A)",
  "lday_leq_db_a",              "$L_\\mathrm{day}$",    "dB(A)",
  "lnight_leq_db_a",            "$L_\\mathrm{night}$",  "dB(A)",
  "lpeak_leq_db_a",             "$L_\\mathrm{peak}$",   "dB(A)",
  "linter_peak_leq_db_a",       "$L_\\mathrm{inter\\text{-}peak}$","dB(A)",
  "mp_leq_db_a",                "Morning peak",         "dB(A)",
  "ep_leq_db_a",                "Evening peak",         "dB(A)",
  "eve_op_leq_db_a",            "Evening off-peak",     "dB(A)",
  "laeq_sd_24h_db_a",           "24 h SD",              "dB",
  "laeq_peak_trough_diff_db",   "Diurnal amplitude",    "dB",
  "high_noise_share_gt_60db",   "Hours $>$60 dB(A)",    "\\%",
  "high_noise_share_gt_65db",   "Hours $>$65 dB(A)",    "\\%")
rows <- metrics %>% rowwise() %>% mutate(
  x = list(nm[[col]]),
  is_share = unit == "\\%",
  mean = ifelse(is_share, f1(mean(x) * 100), f1(mean(x))),
  sd   = ifelse(is_share, f1(sd(x) * 100),   f1(sd(x))),
  med  = ifelse(is_share, f1(median(x) * 100), f1(median(x))),
  mn   = ifelse(is_share, f1(min(x) * 100),   f1(min(x))),
  mx   = ifelse(is_share, f1(max(x) * 100),   f1(max(x)))) %>% ungroup()
body <- rows %>% mutate(line = sprintf("%s (%s) & %s & %s & %s & %s & %s \\\\",
                        lab, unit, mean, sd, med, mn, mx)) %>% pull(line)
wr(c("\\begin{tabular}{l r r r r r}", "\\toprule",
     "Metric & Mean & SD & Median & Min & Max \\\\", "\\midrule",
     body, "\\bottomrule", "\\end{tabular}"), "si_metrics.tex")

# -----------------------------------------------------------------------------
# Sound-event classifier validation against the manually labelled subsample
# -----------------------------------------------------------------------------
if (file.exists(file.path(EXPLORE_OUT, "33_sound_event_validation__summary.csv"))) {
  vs <- read_csv(file.path(EXPLORE_OUT, "33_sound_event_validation__summary.csv"), show_col_types = FALSE)
  vp <- read_csv(file.path(EXPLORE_OUT, "33_sound_event_validation__perclass.csv"), show_col_types = FALSE)
  va <- read_csv(file.path(EXPLORE_OUT, "33_sound_event_validation__shareagreement.csv"), show_col_types = FALSE)
  glab <- c(traffic = "Traffic", human = "Human", natural = "Natural", mechanical = "Mechanical")
  gord <- c("traffic", "human", "natural", "mechanical")
  s1m <- function(x) ifelse(x >= 0, paste0("+", f1(x)), paste0("-", f1(abs(x))))
  arow <- function(res, disp, lv, sb) {
    r <- vs[vs$task == "pri" & vs$level == lv & vs$subset == sb, ]
    sprintf("%s & %s & %d & %s & %s \\\\", res, disp, r$N_valid_majority[1], f3(r$accuracy[1]), f3(r$kappa[1]))
  }
  abody <- c(arow("Four groups", "All", "group", "all"),
             arow("Four groups", "Night", "group", "night"),
             arow("Twelve categories", "All", "code", "all"),
             arow("Twelve categories", "Night", "code", "night"))
  pbody <- vapply(gord, function(g) { r <- vp[vp$group == g, ]
    sprintf("%s & %s & %s & %s & %d \\\\", glab[[g]], f3(r$precision[1]), f3(r$recall[1]), f3(r$f1[1]), r$support_majority[1]) }, character(1))
  sbody <- vapply(gord, function(g) { r <- va[va$group == g, ]
    sprintf("%s & %s & %s & $%s$ & $[%s, %s]$ \\\\", glab[[g]], f1(r$model_share_pct[1]), f1(r$human_majority_share_pct[1]), s1m(r$diff_pp[1]), s1m(r$ci_lo_pp[1]), s1m(r$ci_hi_pp[1])) }, character(1))
  wr(c(
    "\\emph{(a) Agreement with the human reference (dominant source)}\\par\\smallskip",
    "\\begin{tabular}{l l r r r}", "\\toprule",
    "Resolution & Subset & $N$ & Accuracy & Cohen's $\\kappa$ \\\\", "\\midrule", abody, "\\bottomrule",
    "\\end{tabular}\\par\\medskip",
    "\\emph{(b) Per-group reliability (four groups, dominant source, all segments)}\\par\\smallskip",
    "\\begin{tabular}{l r r r r}", "\\toprule",
    "Source group & Precision & Recall & $F_1$ & Segments \\\\", "\\midrule", pbody, "\\bottomrule",
    "\\end{tabular}\\par\\medskip",
    "\\emph{(c) Aggregate composition: model vs human majority (four groups)}\\par\\smallskip",
    "\\begin{tabular}{l r r r r}", "\\toprule",
    "Source group & Model (\\%) & Human (\\%) & Difference (pp) & 95\\% CI (pp) \\\\", "\\midrule", sbody, "\\bottomrule",
    "\\end{tabular}"), "si_validation.tex")
} else {
  message("  si_validation.tex skipped: the classifier-validation results are absent (data/validation/ not supplied)")
}

# -----------------------------------------------------------------------------
# PCA variance + loadings (PC1-PC3)
# -----------------------------------------------------------------------------
load <- read_csv(file.path(EXPLORE_OUT, "18_outcome_structure__pca_loadings.csv"), show_col_types = FALSE)
vexp <- read_csv(file.path(EXPLORE_OUT, "18_outcome_structure__pca_variance.csv"), show_col_types = FALSE)
lvl <- c("LAeq_24h","Lday","Lnight","Lpeak","Linter-peak","Morning peak",
         "Evening peak","Evening off-peak","24h SD","Peak-trough",
         "Share >60 dB","Share >65 dB")
lab_tex <- c("LAeq_24h"="$L_\\mathrm{Aeq,24h}$","Lday"="$L_\\mathrm{day}$",
  "Lnight"="$L_\\mathrm{night}$","Lpeak"="$L_\\mathrm{peak}$",
  "Linter-peak"="$L_\\mathrm{inter\\text{-}peak}$","Morning peak"="Morning peak",
  "Evening peak"="Evening peak","Evening off-peak"="Evening off-peak",
  "24h SD"="24 h SD","Peak-trough"="Diurnal amplitude",
  "Share >60 dB"="Hours $>$60 dB(A)","Share >65 dB"="Hours $>$65 dB(A)")
lo <- load %>% mutate(metric_label = factor(metric_label, levels = lvl)) %>% arrange(metric_label)
body <- lo %>% mutate(line = sprintf("%s & %s & %s & %s \\\\", lab_tex[as.character(metric_label)],
                      sgn(PC1), sgn(PC2), sgn(PC3))) %>% pull(line)
vline <- sprintf("Variance explained & %s\\%% & %s\\%% & %s\\%% \\\\",
                 f1(vexp$variance_explained[1]*100), f1(vexp$variance_explained[2]*100),
                 f1(vexp$variance_explained[3]*100))
wr(c("\\begin{tabular}{l r r r}", "\\toprule",
     "Metric & PC1 & PC2 & PC3 \\\\", "\\midrule", body, "\\midrule", vline,
     "\\bottomrule", "\\end{tabular}"), "si_pca.tex")

# -----------------------------------------------------------------------------
# Sound-source composition correlations (full)
# -----------------------------------------------------------------------------
sc <- read_csv(file.path(EXPLORE_OUT, "08_sound_source__corr_spearman.csv"), show_col_types = FALSE)
o5 <- c(laeq_24h_leq_db_a="$L_\\mathrm{Aeq,24h}$", lday_leq_db_a="$L_\\mathrm{day}$",
  lnight_leq_db_a="$L_\\mathrm{night}$", lpeak_leq_db_a="$L_\\mathrm{peak}$",
  linter_peak_leq_db_a="$L_\\mathrm{inter\\text{-}peak}$",
  laeq_sd_24h_db_a="24 h SD", laeq_peak_trough_diff_db="Diurnal amplitude")
src <- c(sound_traffic_share="Traffic", sound_human_share="Human",
         sound_mechanical_share="Mechanical", sound_natural_share="Natural")
wide <- sc %>% filter(outcome %in% names(o5), predictor %in% names(src)) %>%
  mutate(cell = paste0(sgn(spearman_rho), star(p_value))) %>%
  select(outcome, predictor, cell) %>%
  pivot_wider(names_from = predictor, values_from = cell)
wide <- wide[match(names(o5), wide$outcome), ]
body <- apply(wide, 1, function(r) sprintf("%s & %s & %s & %s & %s \\\\",
        o5[r["outcome"]], r["sound_traffic_share"], r["sound_human_share"],
        r["sound_mechanical_share"], r["sound_natural_share"]))
wr(c("\\begin{tabular}{l r r r r}", "\\toprule",
     "Outcome & Traffic & Human & Mechanical & Natural \\\\", "\\midrule",
     body, "\\bottomrule", "\\end{tabular}"), "si_source_corr.tex")

# -----------------------------------------------------------------------------
# Domain cross-validation performance (all models x outcomes), with the
# paired count of repeats in which the source-only model beat each rival
# -----------------------------------------------------------------------------
cv <- read_csv(file.path(EXPLORE_OUT, "21_domain_cv__model_performance.csv"), show_col_types = FALSE)
gn <- read_csv(file.path(EXPLORE_OUT, "21_domain_cv__model_gain.csv"), show_col_types = FALSE) %>%
  select(outcome, model, cv_rmse_gain_vs_null)
pw <- read_csv(file.path(EXPLORE_OUT, "21_domain_cv__paired_wins.csv"), show_col_types = FALSE) %>%
  select(outcome, model, source_wins)
mod_lab <- c(intercept_only="Intercept only", controls="Controls", environment="Environment",
  activity="Activity", source="Source", environment_source="Environment + source",
  exposure_domains="Activity + environment + source", all_domains_controls="All domains + controls")
o6 <- c(laeq_24h_leq_db_a="$L_\\mathrm{Aeq,24h}$", lday_leq_db_a="$L_\\mathrm{day}$",
  lnight_leq_db_a="$L_\\mathrm{night}$", lpeak_leq_db_a="$L_\\mathrm{peak}$",
  linter_peak_leq_db_a="$L_\\mathrm{inter\\text{-}peak}$")
cj <- cv %>% left_join(gn, by = c("outcome","model")) %>%
  left_join(pw, by = c("outcome","model")) %>%
  filter(outcome %in% names(o6)) %>%
  mutate(outcome = factor(outcome, levels = names(o6)),
         model = factor(model, levels = names(mod_lab))) %>%
  arrange(outcome, model)
bld <- function(s, on) ifelse(on, paste0("\\textbf{", s, "}"), s)
body <- c()
for (oc in names(o6)) {
  body <- c(body, sprintf("\\multicolumn{6}{l}{\\emph{%s}}\\\\", o6[oc]))
  sub <- cj %>% filter(outcome == oc)
  body <- c(body, sub %>% mutate(
    is_src = model == "source",
    wins = ifelse(is_src, "--", as.character(source_wins)),
    line = sprintf("\\quad %s & %s & %s & %s & %s & %s \\\\",
      bld(mod_lab[as.character(model)], is_src), bld(as.character(p), is_src),
      bld(f3(adj_r2), is_src),
      bld(paste0(f2(cv_rmse_mean), "$\\pm$", f2(cv_rmse_sd)), is_src),
      bld(f2(cv_rmse_gain_vs_null), is_src), wins)) %>% pull(line))
}
wr(c("\\begin{tabular}{l r r r r r}", "\\toprule",
     "Outcome / model & Predictors & Adj.\\ $R^2$ & CV RMSE (dB) & \\makecell{$\\Delta$ vs\\\\intercept} & \\makecell{Source lower\\\\(of 100)} \\\\", "\\midrule",
     body, "\\bottomrule", "\\end{tabular}"), "si_domain_cv.tex")

# -----------------------------------------------------------------------------
# Group differences (Kruskal-Wallis) by line and station type
# -----------------------------------------------------------------------------
kw <- read_csv(file.path(EXPLORE_OUT, "04_groupdiff__kruskal_wallis.csv"), show_col_types = FALSE)
mlab <- c(laeq_24h_leq_db_a="$L_\\mathrm{Aeq,24h}$", lday_leq_db_a="$L_\\mathrm{day}$",
  lnight_leq_db_a="$L_\\mathrm{night}$", lpeak_leq_db_a="$L_\\mathrm{peak}$",
  linter_peak_leq_db_a="$L_\\mathrm{inter\\text{-}peak}$",
  laeq_peak_trough_diff_db="Diurnal amplitude", laeq_sd_24h_db_a="24 h SD")
glab2 <- c(line_name = "Metro line (5 groups)", station_type_code = "Station type (3 groups)")
ko <- kw %>% filter(metric %in% names(mlab)) %>%
  mutate(grouping = factor(grouping, levels = c("line_name","station_type_code")),
         metric = factor(metric, levels = names(mlab))) %>%
  arrange(grouping, metric)
body <- c()
for (g in levels(ko$grouping)) {
  body <- c(body, sprintf("\\multicolumn{5}{l}{\\emph{%s}}\\\\", glab2[g]))
  sub <- ko %>% filter(grouping == g)
  body <- c(body, sub %>% mutate(line = sprintf("\\quad %s & %s & %d & %s & %s \\\\",
            mlab[as.character(metric)], f2(h_stat), k_groups - 1, pfmt(p_value), f3(epsilon2))) %>% pull(line))
}
wr(c("\\begin{tabular}{l r r r r}", "\\toprule",
     "Grouping / metric & $H$ & df & $p$ & $\\epsilon^2$ \\\\", "\\midrule",
     body, "\\bottomrule", "\\end{tabular}"), "si_groupdiff.tex")

# -----------------------------------------------------------------------------
# Source-association influence robustness (jackknife / bootstrap / subsets)
# -----------------------------------------------------------------------------
jk <- read_csv(file.path(EXPLORE_OUT, "27_source_robustness__jackknife_summary.csv"), show_col_types = FALSE)
bc <- read_csv(file.path(EXPLORE_OUT, "27_source_robustness__boot_ci.csv"), show_col_types = FALSE)
ss <- read_csv(file.path(EXPLORE_OUT, "27_source_robustness__subsets.csv"), show_col_types = FALSE)
rob <- jk %>% select(label, full_rho, loo_min, loo_max, n_sign_flip, n_lose_sig) %>%
  left_join(bc %>% select(label, ci_lo, ci_hi), by = "label") %>%
  left_join(ss %>% select(label, rho_excl_zero_natural, rho_trim5_traffic), by = "label")
labtex <- c("traffic_share ~ LAeq_24h"="Traffic $\\sim L_\\mathrm{Aeq,24h}$",
  "natural_share ~ LAeq_24h"="Natural $\\sim L_\\mathrm{Aeq,24h}$",
  "traffic_share ~ Lnight"="Traffic $\\sim L_\\mathrm{night}$",
  "natural_share ~ Lnight"="Natural $\\sim L_\\mathrm{night}$")
body <- rob %>% mutate(line = sprintf("%s & %s & [%s, %s] & [%s, %s] & %d / %d & %s & %s \\\\",
        labtex[label], sgn(full_rho), f2(loo_min), f2(loo_max), f2(ci_lo), f2(ci_hi),
        n_sign_flip, n_lose_sig, sgn(rho_excl_zero_natural), sgn(rho_trim5_traffic))) %>% pull(line)
wr(c("\\begin{tabular}{l r c c c r r}", "\\toprule",
     "Association & $\\rho$ & LOO range & Boot.\\ 95\\% CI & \\makecell{Flips /\\\\lose sig.} & \\makecell{Excl.\\\\zero-nat.} & \\makecell{Trim 5\\%\\\\traffic} \\\\", "\\midrule",
     body, "\\bottomrule", "\\end{tabular}"), "si_robustness.tex")

# -----------------------------------------------------------------------------
# Built-environment FDR-significant partial-Spearman associations (longtable rows)
# -----------------------------------------------------------------------------
be <- read_csv(file.path(EXPLORE_OUT, "06_built_env__correlations.csv"), show_col_types = FALSE)
catlab <- c(road_arterial_density="Arterial road density", poi_total_density="POI density",
  commercial_land_pct="Commercial land", landuse_mix_shannon="Land-use mix",
  residential_land_pct="Residential land", road_pedestrian_density="Pedestrian road density",
  road_expressway_density="Expressway density", green_land_pct="Green-space land",
  water_body_pct="Water body")
olab <- c(laeq_24h_leq_db_a="$L_\\mathrm{Aeq,24h}$", lday_leq_db_a="$L_\\mathrm{day}$",
  lnight_leq_db_a="$L_\\mathrm{night}$", lpeak_leq_db_a="$L_\\mathrm{peak}$",
  lpeak_db_a="$L_\\mathrm{peak\\,inst}$", linter_peak_leq_db_a="$L_\\mathrm{inter\\text{-}peak}$",
  loffpeak_day_leq_db_a="$L_\\mathrm{offpeak,day}$")
sig <- be %>% filter(method == "partial_spearman", q_fdr < 0.05) %>%
  arrange(q_fdr) %>%
  mutate(o = olab[outcome], c = catlab[category]) %>%
  filter(!is.na(o), !is.na(c))
body <- sig %>% mutate(line = sprintf("%s & %s & %d m & %s & %s \\\\",
        o, c, scale_m, sgn(rho), f3(q_fdr))) %>% pull(line)
wr(body, "si_builtenv_rows.tex")

# -----------------------------------------------------------------------------
# Spatial source-term coefficients + residual Moran's I
# -----------------------------------------------------------------------------
st <- read_csv(file.path(EXPLORE_OUT, "28_source_spatial__source_term_summary.csv"), show_col_types = FALSE)
mo <- read_csv(file.path(EXPLORE_OUT, "28_source_spatial__residual_moran.csv"), show_col_types = FALSE)
tlab <- c(sound_traffic_share="Traffic share", sound_natural_share="Natural-source share")
oo <- c("LAeq_24h"="$L_\\mathrm{Aeq,24h}$","Lday"="$L_\\mathrm{day}$","Lnight"="$L_\\mathrm{night}$")
stx <- st %>% mutate(outcome = factor(outcome, levels = names(oo)),
                     term = factor(term, levels = names(tlab))) %>% arrange(outcome, term)
cbody <- stx %>% mutate(line = sprintf("%s & %s & %s & %s & %s & %s \\\\",
          oo[as.character(outcome)], tlab[as.character(term)],
          f1(coef__OLS), f1(coef__OLS_trend), f1(coef__SAR), f1(coef__SEM))) %>% pull(line)
spm <- mo %>% mutate(spec = recode(spec,
  "intercept_only"="int","OLS[sound_traffic_share]"="traf","OLS_trend[sound_traffic_share]"="traf_tr",
  "OLS[sound_natural_share]"="nat","OLS_trend[sound_natural_share]"="nat_tr")) %>%
  select(outcome, spec, resid_moran_I) %>% pivot_wider(names_from = spec, values_from = resid_moran_I)
spm <- spm[match(names(oo), spm$outcome), ]
mbody <- apply(spm, 1, function(r) sprintf("%s & %s & %s & %s & %s & %s \\\\",
          oo[r["outcome"]], f2(as.numeric(r["int"])), f2(as.numeric(r["traf"])),
          f2(as.numeric(r["traf_tr"])), f2(as.numeric(r["nat"])), f2(as.numeric(r["nat_tr"]))))
wr(c("\\begin{tabular}{l l r r r r}", "\\toprule",
     "\\multicolumn{6}{l}{\\emph{Source-share coefficient (dB per unit share)}}\\\\",
     "Outcome & Term & OLS & Trend & Spatial lag & Spatial error \\\\", "\\midrule",
     cbody, "\\midrule",
     "\\multicolumn{6}{l}{\\emph{Residual Moran's I}}\\\\",
     "Outcome & Intercept & +Traffic & +Traffic, trend & +Natural & +Natural, trend \\\\", "\\midrule",
     mbody, "\\bottomrule", "\\end{tabular}"), "si_spatial.tex")

# -----------------------------------------------------------------------------
# Beyond-loudness partial correlations (natural share vs NDVI under controls)
# -----------------------------------------------------------------------------
bl <- read_csv(file.path(EXPLORE_OUT, "29_beyond_loudness_partial.csv"), show_col_types = FALSE)
ylab <- c(sound_natural_share="Natural-source share")
clab <- c(none="zero-order", laeq="$\\mid L_\\mathrm{Aeq,24h}$", lnight="$\\mid L_\\mathrm{night}$",
  laeq_plus_other="$\\mid L_\\mathrm{Aeq}$ + road", lnight_plus_other="$\\mid L_\\mathrm{night}$ + road")
blf <- bl %>% filter(subset == "all_109", predictor == "ndvi_500",
                     control %in% names(clab)) %>%
  mutate(y = ylab[outcome], control = factor(control, levels = names(clab))) %>%
  filter(!is.na(y)) %>% arrange(outcome, control) %>%
  mutate(cell = paste0(sgn(rho)))
wide <- blf %>% select(y, outcome, control, cell) %>%
  pivot_wider(names_from = control, values_from = cell)
wide <- wide[match(names(ylab), wide$outcome), ]
body <- apply(wide, 1, function(r) sprintf("%s & %s & %s & %s & %s & %s \\\\",
        r["y"], r["none"], r["laeq"], r["lnight"], r["laeq_plus_other"], r["lnight_plus_other"]))
wr(c("\\begin{tabular}{l r r r r r}", "\\toprule",
     "Measure vs NDVI 500\\,m & Zero-order & $\\mid L_\\mathrm{Aeq}$ & $\\mid L_\\mathrm{night}$ & $\\mid L_\\mathrm{Aeq}$+road & $\\mid L_\\mathrm{night}$+road \\\\", "\\midrule",
     body, "\\bottomrule", "\\end{tabular}"), "si_beyond_loudness.tex")

message("All SI table fragments generated in ", TBL_DIR)
