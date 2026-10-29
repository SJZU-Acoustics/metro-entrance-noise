suppressPackageStartupMessages(library(tidyverse))
EXPLORE_OUT <- "intermediate/analysis"
TBL_DIR <- "intermediate/support_tables"
dir.create(TBL_DIR,recursive=TRUE,showWarnings=FALSE)
mm <- function(s) sub("^-", "$-$", s)
f1 <- function(x) mm(formatC(x, format = "f", digits = 1))
f2 <- function(x) mm(formatC(x, format = "f", digits = 2))
f3 <- function(x) mm(formatC(x, format = "f", digits = 3))
sgn <- function(x) ifelse(x >= 0, paste0("$+$", f2(x)), paste0("$-$", f2(abs(x))))
pfmt <- function(p) ifelse(p < 0.001, "$<$0.001", formatC(p, format = "f", digits = 3))
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
wr <- function(lines, file) { writeLines(lines, file.path(TBL_DIR, file));
  message("  wrote ", file) }

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
