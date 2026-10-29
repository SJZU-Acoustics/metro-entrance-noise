suppressPackageStartupMessages(library(tidyverse))
script <- "code/build_displays.R"
w <- "output"; project <- "."; e <- "intermediate"
source('code/shared/helpers.R'); source('code/shared/plot_style.R')
d <- load_p11(e)
a37 <- 'analysis_37_information_choices'; a38 <- 'analysis_38_priority_sources'; a39 <- 'analysis_39_persistent_error'
a40 <- 'analysis_40_assessment_benchmarks'; a41 <- 'analysis_41_priority_profile_reference'
a42 <- 'analysis_42_priority_places'; a43 <- 'analysis_43_recording_dates'; a35 <- 'analysis_35_assessment_information'
used <- character()
rd <- function(a, f) { p <- file.path(e, a, 'outcome', paste0(f, '.csv')); used <<- union(used, p); read_csv(p, show_col_types = FALSE) }
fig <- file.path(w, 'figures'); tab <- file.path(w, 'tables'); out <- 'intermediate/display_data'
blue <- '#0072B2'; orange <- '#E69F00'; green <- '#009E73'; red <- '#D55E00'; sky <- '#56B4E9'; pink <- '#CC79A7'
three <- theme(axis.title = element_text(size = 9), axis.text = element_text(size = 8), legend.text = element_text(size = 8))

# Figure 1. Screening procedure; blocks follow the three Results subsections.
box <- tibble(xmin = c(4, 4, 62, 120, 4), xmax = c(174, 58, 116, 174, 174), ymin = c(64, 16, 16, 16, 1), ymax = c(77, 56, 56, 56, 10))
steps <- tibble(x = c(31, 89, 147), head = c('RQ1  Select entrances', 'RQ2  Characterise sources', 'RQ3  Qualify the list'),
  ask = c('Which entrances, in which period?', 'Which sources, at which hours?', 'How stable is the selection?'),
  items = c('24-hour and period indices\nLevel, activity and joint rules\nHourly or period alignment\nList sizes 1–30\nFixed limit values',
    'Day and night priority lists\nHourly source groups\nWhole network and random sets\nList sizes 1–30',
    'Errors shared across hours\nErrors shared across dates\nOmitted entrances and dates'))
arrow_small <- grid::arrow(length = grid::unit(1.6, 'mm'), type = 'closed')
g <- ggplot() + geom_rect(data = box[1:4, ], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white', colour = 'black', linewidth = .35) +
  geom_rect(data = box[5, ], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = 'white', colour = 'black', linewidth = .35, linetype = 2) +
  annotate('text', x = 89, y = 73, label = 'Network survey', fontface = 'bold', size = 9 / .pt, family = 'Helvetica') +
  annotate('text', x = 89, y = 67.6, label = '109 metro entrances, 24 hourly records each: sound level, activity proxy and classified sound source', size = 8.5 / .pt, family = 'Helvetica') +
  geom_text(data = steps, aes(x, 51.5, label = head), fontface = 'bold', size = 9 / .pt, family = 'Helvetica') +
  geom_text(data = steps, aes(x, 45.5, label = ask), fontface = 'italic', size = 8.5 / .pt, family = 'Helvetica') +
  geom_text(data = steps, aes(x - 23, 41, label = items), hjust = 0, vjust = 1, lineheight = 1.12, size = 8.5 / .pt, family = 'Helvetica') +
  geom_segment(data = steps, aes(x = x, xend = x, y = 64, yend = 56.4), arrow = arrow_small, linewidth = .35) +
  annotate('segment', x = c(58, 116), xend = c(61.6, 119.6), y = 36, yend = 36, arrow = arrow_small, linewidth = .35) +
  geom_segment(data = steps, aes(x = x, xend = x, y = 16, yend = 10.4), arrow = arrow_small, linewidth = .35) +
  annotate('text', x = 89, y = 5.5, label = 'Screening list with a source profile and a stated stability, for site investigation', size = 8.5 / .pt, family = 'Helvetica') +
  coord_cartesian(xlim = c(0, 178), ylim = c(0, 78), expand = FALSE) + theme_void() + theme(plot.background = element_rect(fill = 'white', colour = NA))
ggsave(file.path(fig, 'framework.png'), g, width = 178 / 25.4, height = 78 / 25.4, units = 'in', dpi = 600, bg = 'white')

# Figure 2. Existing road cache is immutable background; data remain canonical.
sm <- read_csv(file.path(e, 'data/clean/station_master.csv'), show_col_types = FALSE)
map <- sm %>% transmute(station_key, x = lng_round, y = lat_round, valid = station_key %in% d$master$station_key)
roadpath <- 'data/osm_roads_shenyang.csv'
roads <- read_csv(roadpath, show_col_types = FALSE)
xr <- range(map$x) + c(-.01, .01); yr <- range(map$y) + c(-.01, .01)
span <- diff(yr) / cos(mean(map$y) * pi / 180); xr <- mean(xr) + c(-.5, .5) * span
roads <- filter(roads, lng > xr[1], lng < xr[2], lat > yr[1], lat < yr[2])
base_map <- function() ggplot(map, aes(x, y)) + geom_path(data = roads, aes(lng, lat, group = line_id), inherit.aes = FALSE, colour = 'grey85', linewidth = .2) +
  coord_fixed(ratio = 1 / cos(mean(map$y) * pi / 180), xlim = xr, ylim = yr) + scale_x_continuous(breaks = c(123.2, 123.4, 123.6)) +
  labs(x = 'Longitude (°E)', y = 'Latitude (°N)') + p11_theme() + theme(aspect.ratio = NULL)
pa <- base_map() + geom_point(aes(shape = valid), size = 1.2, colour = blue) +
  scale_shape_manual(values = c('FALSE' = 1, 'TRUE' = 16), labels = c('No valid record', 'Valid record')) + theme(legend.key.width = grid::unit(3, 'mm'), legend.key.spacing.x = grid::unit(1, 'mm'))
hr <- d$panel %>% group_by(hour) %>% summarise(mean = mean(laeq_hour_db_a), q25 = quantile(laeq_hour_db_a, .25), q75 = quantile(laeq_hour_db_a, .75), .groups = 'drop')
pb <- ggplot(hr, aes(hour, mean)) + geom_ribbon(aes(ymin = q25, ymax = q75), fill = blue, alpha = .18) + geom_line(colour = blue, linewidth = .6) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 23)) + labs(x = 'Hour of day', y = 'Hourly sound level (dB(A))') + p11_theme()
ind <- rd(a40, 'indicators')
lev <- ind %>% select(station_key, Day = lday, Night = lnight) %>% pivot_longer(-station_key, names_to = 'period', values_to = 'level')
pc <- ggplot(lev, aes(level, colour = period, linetype = period)) + geom_vline(xintercept = c(55, 70), colour = c(blue, orange), linewidth = .35, linetype = 3) +
  stat_bin(geom = 'step', binwidth = 2, boundary = 0, position = 'identity', linewidth = .6, direction = 'mid') +
  scale_colour_manual(values = c(Day = orange, Night = blue)) + scale_linetype_manual(values = c(Day = 1, Night = 2)) +
  scale_x_continuous(limits = c(44, 80), breaks = seq(45, 75, 10)) + scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = 'Period sound level (dB(A))', y = 'Number of entrances') + p11_theme() + theme(legend.key.width = grid::unit(7, 'mm'))
p11_png(((pa | pb | pc) + plot_annotation(tag_levels = 'a')) & three, file.path(fig, 'study.png'), 76)

# Figure 3 and Figure S1. Membership and coverage against the period-joint list at every size.
x <- rd(a37, 'coverage')
ref <- x %>% filter(method == 'period_joint') %>% select(definition, period, budget, reference = fraction_index)
x <- left_join(x, ref, by = c('definition', 'period', 'budget')) %>% mutate(gap = 100 * (reference - fraction_index), shared = overlap / budget)
cols <- c(full_level = blue, period_level = sky, full_joint = red, period_product = green)
labs_m <- c(full_level = '24 h level', period_level = 'Period level', full_joint = '24 h joint', period_product = 'Period product')
plab <- c(day = 'Day', night = 'Night', peak = 'Peak', offpeak = 'Inter-peak')
curve <- function(p, metric) ggplot(filter(x, definition == 'raw50', period == p, method %in% names(cols)) %>% mutate(method = factor(method, names(cols))),
    aes(budget, .data[[metric]], colour = method, linetype = method)) +
  geom_vline(xintercept = 10, colour = 'grey70', linewidth = .25) + geom_line(linewidth = .6) +
  scale_colour_manual(values = cols, labels = labs_m) + scale_linetype_manual(values = c(full_level = 2, period_level = 3, full_joint = 1, period_product = 4), labels = labs_m) +
  scale_x_continuous(breaks = c(1, 5, 10, 20, 30)) +
  labs(x = 'Number of selected entrances', y = paste0(plab[[p]], ifelse(metric == 'gap', ': index coverage gap (pp)', ': shared fraction of list'))) +
  (if (metric == 'gap') scale_y_continuous(limits = c(0, 15), breaks = c(0, 3, 6, 9, 12, 15), expand = expansion(mult = c(0, .03))) else scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25))) + p11_theme()
four <- function(p1, p2) ((curve(p1, 'shared') | curve(p2, 'shared')) / (curve(p1, 'gap') | curve(p2, 'gap'))) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a')
p11_png(four('day', 'night') & theme(legend.position = 'top'), file.path(fig, 'selection.png'))
p11_png(four('peak', 'offpeak') & theme(legend.position = 'top'), file.path(fig, 'si/other_periods.png'))

# Figure 4. Shared, day-only and night-only priorities on the map and against their day and night levels.
pl <- rd(a42, 'priority_places') %>% mutate(class = factor(priority_class, c('other', 'shared', 'day_only', 'night_only')))
ccol <- c(other = 'grey65', shared = 'black', day_only = orange, night_only = blue); cshape <- c(other = 16, shared = 16, day_only = 17, night_only = 15)
csize <- c(other = .8, shared = 2, day_only = 2, night_only = 2); clab <- c(other = 'Other', shared = 'Day and night', day_only = 'Day only', night_only = 'Night only')
class_scales <- list(scale_colour_manual(values = ccol, labels = clab), scale_shape_manual(values = cshape, labels = clab), scale_size_manual(values = csize, labels = clab))
pl <- arrange(pl, class)
m1 <- base_map() + geom_point(data = pl, aes(lng_round, lat_round, colour = class, shape = class, size = class)) + class_scales
m2 <- ggplot(pl, aes(lday, lnight, colour = class, shape = class, size = class)) + geom_abline(slope = 1, intercept = 0, colour = 'grey50', linewidth = .3) +
  geom_point() + class_scales + scale_x_continuous(limits = c(46, 77), breaks = seq(50, 75, 5)) + scale_y_continuous(limits = c(46, 77), breaks = seq(50, 75, 5)) +
  labs(x = 'Daytime sound level (dB(A))', y = 'Night-time sound level (dB(A))') + p11_theme()
p11_png(((m1 | m2) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a')) & theme(legend.position = 'top'), file.path(fig, 'priority_places.png'), 105)

# Figure 5. Hourly source tiles, as produced and inspected in analysis 38.
tiles <- file.path(e, a38, 'outcome/priority_source_hours.png'); used <- union(used, tiles)
stopifnot(file.copy(tiles, file.path(fig, 'source_hours.png'), overwrite = TRUE))

# Figure 6. Own-period traffic share of the priority lists against the network and random sets of the same size.
pr <- rd(a41, 'profile_by_size') %>% filter(source == 'traffic', rule %in% c('joint', 'level'))
rr <- rd(a41, 'random_reference') %>% filter(source == 'traffic'); net <- rd(a41, 'network_reference') %>% filter(source == 'traffic')
prof_panel <- function(p) {
  z <- pr %>% transmute(rule = factor(rule, c('joint', 'level')), budget, share = 100 * if (p == 'day') day_day else night_night)
  ggplot(z, aes(budget, share)) + geom_ribbon(data = filter(rr, period == p), aes(budget, ymin = 100 * random_q025, ymax = 100 * random_q975), inherit.aes = FALSE, fill = 'grey85') +
    geom_hline(yintercept = 100 * net$share[net$period == p], linewidth = .35, linetype = 3) + geom_vline(xintercept = 10, colour = 'grey70', linewidth = .25) +
    geom_line(aes(colour = rule, linetype = rule), linewidth = .6) + scale_colour_manual(values = c(joint = red, level = blue), labels = c('Joint list', 'Level list')) +
    scale_linetype_manual(values = c(joint = 1, level = 2), labels = c('Joint list', 'Level list')) + scale_x_continuous(breaks = c(1, 5, 10, 20, 30)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, .02))) +
    labs(x = 'Number of selected entrances', y = paste0(plab[[p]], ': traffic share of hours (%)')) + p11_theme()
}
p11_png(((prof_panel('day') | prof_panel('night')) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a')) & theme(legend.position = 'top'), file.path(fig, 'source_profile.png'), 100)

# Figure 7. Retained priorities under persistent error, and how often each entrance is selected.
sim <- rd(a39, 'simulation_summary'); inc <- rd(a39, 'inclusion_frequencies') %>% filter(scenario == 'main_rho1')
elab <- c('1 dB / 20%', '2 dB / 30%')
err_panel <- function(p) ggplot(filter(sim, period == p), aes(persistence, observed_overlap_mean, colour = scale, shape = scale)) +
  geom_errorbar(aes(ymin = observed_overlap_q05, ymax = observed_overlap_q95), width = .04, position = position_dodge(.08), linewidth = .4) +
  geom_line(aes(group = scale), position = position_dodge(.08), linewidth = .5) + geom_point(position = position_dodge(.08), size = 2) +
  scale_colour_manual(values = c(main = blue, stress = red), labels = elab) + scale_shape_manual(values = c(main = 16, stress = 17), labels = elab) +
  scale_x_continuous(breaks = c(0, .5, 1), labels = c('0', '50', '100'), limits = c(-.06, 1.06)) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  labs(x = 'Error variance shared across hours (%)', y = paste0(plab[[p]], ': original priorities retained (of 10)')) + p11_theme()
inc_panel <- function(p) {
  z <- inc %>% filter(period == p, inclusion_frequency >= .05) %>% left_join(select(pl, station_key, rank_day, rank_night), by = 'station_key') %>%
    mutate(listed = if (p == 'day') rank_day <= 10 else rank_night <= 10) %>% arrange(inclusion_frequency) %>% mutate(station_key = factor(station_key, station_key))
  ggplot(z, aes(100 * inclusion_frequency, station_key)) + geom_segment(aes(x = 0, xend = 100 * inclusion_frequency, yend = station_key), colour = 'grey70', linewidth = .3) +
    geom_point(aes(shape = listed), colour = blue, size = 1.8) + scale_shape_manual(values = c('TRUE' = 16, 'FALSE' = 1), labels = c('TRUE' = 'Original priority', 'FALSE' = 'Other entrance'), breaks = c('TRUE', 'FALSE')) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, .03))) +
    labs(x = paste0(plab[[p]], ': selected in simulations (%)'), y = 'Entrance') + p11_theme() + theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size = 8))
}
g <- ((err_panel('day') | err_panel('night')) + plot_layout(guides = 'collect')) / ((inc_panel('day') | inc_panel('night')) + plot_layout(guides = 'collect')) + plot_annotation(tag_levels = 'a')
p11_png(g & theme(legend.position = 'top', legend.text = element_text(size = 8)), file.path(fig, 'reliability.png'), 200)

# Figure S2. Period-product coverage gap under the four activity definitions.
dcols <- c(raw50 = blue, raw20 = red, sqrt50 = green, cap50 = pink); dlab <- c(raw50 = '50 m', raw20 = '20 m', sqrt50 = 'Square-root', cap50 = 'Capped')
prod_panel <- function(p) ggplot(filter(x, period == p, method == 'period_product'), aes(budget, gap, colour = definition, linetype = definition)) +
  geom_vline(xintercept = 10, colour = 'grey70', linewidth = .25) + geom_line(linewidth = .6) + scale_colour_manual(values = dcols, labels = dlab) +
  scale_linetype_manual(values = c(1, 2, 3, 4), labels = dlab) + scale_x_continuous(breaks = c(1, 5, 10, 20, 30)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, .5, 1, 1.5), expand = expansion(mult = c(0, .03))) +
  labs(x = 'Number of selected entrances', y = paste0(plab[[p]], ': product coverage gap (pp)')) + p11_theme()
p11_png(((prod_panel('day') | prod_panel('night')) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a')) & theme(legend.position = 'top'), file.path(fig, 'si/product_definitions.png'), 100)

# Plain table bodies; captions stay in the display scaffolds.
tex <- function(dat, path, align = NULL) {
  if (is.null(align)) align <- paste0('l', paste(rep('r', ncol(dat) - 1), collapse = ''))
  dat[] <- lapply(dat, function(v) { v <- as.character(v); neg <- grepl('^-[0-9]+(\\.[0-9]+)?$', v); v[neg] <- paste0('$-$', sub('^-', '', v[neg])); v })
  lines <- c(paste0('\\begin{tabular}{', align, '}'), '\\toprule', paste0(paste(names(dat), collapse = ' & '), ' \\\\'), '\\midrule',
    apply(dat, 1, function(r) paste0(paste(r, collapse = ' & '), ' \\\\')), '\\bottomrule', '\\end{tabular}')
  writeLines(lines, path)
}
fmt <- function(v, k = 2) { r <- sign(v) * floor(abs(v) * 10^k + .5 + 1e-8) / 10^k; r[r == 0] <- 0; sprintf(paste0('%.', k, 'f'), r) }
rng <- function(a, b, k = 1) ifelse(fmt(a, k) == fmt(b, k), fmt(a, k), paste0(fmt(a, k), '--', fmt(b, k)))
methods <- c('full_level', 'period_level', 'lden_level', 'period_activity', 'full_joint', 'day_joint', 'period_product', 'period_joint')
methodlabel <- c('24 h level', 'Period level', '$L_\\mathrm{den}$', 'Period activity', '24 h joint', 'Daytime joint', 'Period product', 'Period joint')
rules <- rd(a40, 'rules')
rule_table <- function(def, ms = methods, labs = methodlabel) { v <- filter(rules, definition == def, budget == 10)
  as.data.frame(do.call(rbind, lapply(seq_along(ms), function(j) { a <- filter(v, method == ms[j], period == 'day'); b <- filter(v, method == ms[j], period == 'night')
    c(labs[j], a$overlap, fmt(100 * a$fraction_index), b$overlap, fmt(100 * b$fraction_index)) }))) %>%
    setNames(c('Selection rule', 'Day shared', 'Day coverage (\\%)', 'Night shared', 'Night coverage (\\%)')) }
tex(rule_table('primary'), file.path(tab, 'selection.tex'))
tex(rule_table('statutory', c(methods, 'other_definition_joint'), c(methodlabel, 'Primary-period joint')), file.path(tab, 'si/period_definition.tex'))
ex <- rd(a40, 'exceedance') %>% select(reference, limit_period, value_db, definition, n_above) %>% pivot_wider(names_from = definition, values_from = n_above)
tex(data.frame(Reference = ex$reference, Indicator = recode(ex$limit_period, day = 'Day level', night = 'Night level', lden = '$L_\\mathrm{den}$'), Value = ex$value_db,
  Primary = ex$primary, Statutory = ifelse(is.na(ex$statutory), '--', ex$statutory)) %>% setNames(c('Reference', 'Indicator', 'Value (dB)', 'Above, primary', 'Above, statutory')), file.path(tab, 'si/benchmarks.tex'), 'llrrr')
src_order <- c('traffic', 'human', 'mechanical', 'natural')
de <- rd(a41, 'profile_by_size') %>% filter(rule == 'joint', budget == 10) %>% arrange(match(source, src_order))
tex(data.frame(Source = str_to_title(de$source), `Network day` = fmt(100 * de$network_day, 1), `Network night` = fmt(100 * de$network_night, 1), DD = fmt(100 * de$day_day, 1), DN = fmt(100 * de$day_night, 1),
  ND = fmt(100 * de$night_day, 1), NN = fmt(100 * de$night_night, 1), Time = fmt(100 * de$temporal_component, 1), Membership = fmt(100 * de$membership_component, 1),
  Total = fmt(100 * de$total_difference, 1), check.names = FALSE), file.path(tab, 'decomposition.tex'))
ps <- rd(a41, 'profile_by_size') %>% filter(rule == 'joint', source == 'traffic', budget %in% c(5, 10, 15, 20, 30))
rq <- rd(a41, 'random_reference') %>% filter(source == 'traffic', budget %in% ps$budget) %>% transmute(budget, period, r = rng(100 * random_q025, 100 * random_q975)) %>% pivot_wider(names_from = period, values_from = r)
tex(data.frame(Size = ps$budget, DD = fmt(100 * ps$day_day, 1), DN = fmt(100 * ps$day_night, 1), ND = fmt(100 * ps$night_day, 1), NN = fmt(100 * ps$night_night, 1), Time = fmt(100 * ps$temporal_component, 1),
  Membership = fmt(100 * ps$membership_component, 1), Total = fmt(100 * ps$total_difference, 1), `Random, day` = rq$day[match(ps$budget, rq$budget)], `Random, night` = rq$night[match(ps$budget, rq$budget)], check.names = FALSE), file.path(tab, 'si/profile_size.tex'))
# Period product at ten entrances under each activity definition, with the entrance-level alignment ratio.
al <- rd(a40, 'alignment_summary')
a <- filter(x, budget == 10, method == 'period_product') %>% left_join(al, by = c('definition', 'period')) %>% arrange(match(definition, c('raw50', 'raw20', 'sqrt50', 'cap50')), match(period, c('day', 'peak', 'offpeak', 'night')))
tex(data.frame(Definition = unname(dlab[a$definition]), Period = unname(plab[a$period]), Shared = a$overlap, Joint = fmt(100 * a$reference), Product = fmt(100 * a$fraction_index), Gap = fmt(a$gap, 3),
  Ratio = fmt(a$ratio_median), Range = rng(a$ratio_min, a$ratio_max, 2)), file.path(tab, 'si/product.tex'), 'llrrrrrr')
cs <- rd(a42, 'class_summary')
vars <- c(lday = 'Daytime level (dB(A))', lnight = 'Night-time level (dB(A))', night_minus_day = 'Night minus day (dB)', activity_day = 'Daytime activity index ($10^3$)', activity_night = 'Night-time activity index ($10^3$)',
  human_day = 'Human share of daytime hours (\\%)', traffic_night = 'Traffic share of night hours (\\%)', `100m_urban_arterial_road_density_m_m` = 'Arterial-road density (km/km$^2$)',
  `500m_commercial_services_land_pct` = 'Commercial land, 500 m (\\%)', `500m_residential_land_pct` = 'Residential land, 500 m (\\%)', poi_total_density_300m = 'POI density index, 300 m ($10^3$ units)', `500_ndvi` = 'NDVI, 500 m')
scl <- c(activity_day = .001, activity_night = .001, human_day = 100, traffic_night = 100, poi_total_density_300m = .001); dec <- c(`500_ndvi` = 2, poi_total_density_300m = 2)
cell <- function(v, cl) { z <- filter(cs, variable == v, priority_class == cl); s <- if (v %in% names(scl)) scl[[v]] else 1; k <- if (v %in% names(dec)) dec[[v]] else 1
  fm <- function(u) sub('^-', '$-$', fmt(s * u, k)); paste0(fm(z$median), ' (', fm(z$minimum), ' to ', fm(z$maximum), ')') }
ct <- data.frame(Variable = unname(vars), check.names = FALSE)
for (cl in c('shared', 'day_only', 'night_only', 'other')) ct[[c(shared = 'Day and night', day_only = 'Day only', night_only = 'Night only', other = 'Other')[[cl]]]] <- vapply(names(vars), cell, character(1), cl = cl)
writeLines(c('\\begin{tabular}{lrrrr}', '\\toprule', paste0(paste(names(ct), collapse = ' & '), ' \\\\'), '\\midrule', apply(ct, 1, function(r) paste0(paste(r, collapse = ' & '), ' \\\\')), '\\bottomrule', '\\end{tabular}'), file.path(tab, 'si/priority_classes.tex'))
m <- rd(a35, 'disjoint_summary') %>% filter(comparison %in% c('L24_1dB', 'day_night_1dB', 'L24_context', 'day_night_context'))
ml <- c(L24_1dB = '24 h level', day_night_1dB = '24 h, day and night levels', L24_context = '24 h level and context', day_night_context = 'All three levels and context')
tex(data.frame(Matching = unname(ml[m$comparison]), Pairs = m$pairs, Median = fmt(m$median_tv_hours, 1), Six = m$tv_ge6), file.path(tab, 'si/matching.tex'))
su <- rd(a38, 'source_summary') %>% filter(sample == 'all109', label == 'source_group', period %in% c('day', 'night')) %>% arrange(period, match(source, src_order))
tex(data.frame(Period = str_to_title(su$period), Source = str_to_title(su$source), Hours = su$total_hours, Sites = su$sites_present, Share = fmt(100 * su$mean_share, 1), LOO = rng(100 * su$loo_min, 100 * su$loo_max)), file.path(tab, 'si/source_sensitivity.tex'), 'llrrrr')
# Error scenarios: persistent component shared within an entrance (analysis 39) or within a recording date (analysis 43).
s1 <- sim %>% filter(period %in% c('day', 'night')) %>% transmute(within = 'Entrance', scale, share = persistence, period, observed_overlap_mean, observed_overlap_q05, observed_overlap_q95, observed_list_gap_pp_mean, product_gap_pp_mean)
s2 <- rd(a43, 'date_error_summary') %>% transmute(within = 'Date', scale, share = date_share, period, observed_overlap_mean, observed_overlap_q05, observed_overlap_q95, observed_list_gap_pp_mean, product_gap_pp_mean)
s <- bind_rows(s1, s2) %>% arrange(match(within, c('Entrance', 'Date')), scale, share, period)
tex(data.frame(Within = s$within, Error = ifelse(s$scale == 'main', '1 dB / 20\\%', '2 dB / 30\\%'), Shared = sprintf('%.0f\\%%', 100 * s$share), Period = str_to_title(s$period), Retained = fmt(s$observed_overlap_mean),
  Range = paste0(s$observed_overlap_q05, '--', s$observed_overlap_q95), Loss = fmt(s$observed_list_gap_pp_mean), Product = fmt(s$product_gap_pp_mean)) %>% setNames(c('Shared within', 'Error', 'Shared', 'Period', 'Retained', 'Range', 'Loss', 'Product gap')), file.path(tab, 'si/errors.tex'), 'lllrrrrr')
# Omission checks: one listed entrance or one recording date removed, then reselection.
fs <- rd(a43, 'full_sample'); oe <- rd(a43, 'omit_one_entrance'); od <- rd(a43, 'omit_one_date')
q <- tribble(~col, ~label, ~k, 'night_gain_over_full', 'Night gain over 24 h joint list (pp)', 2, 'night_shared_with_full', 'Night entrances shared with 24 h joint list', 0,
  'night_shared_with_level', 'Night entrances shared with night level list', 0, 'night_product_gap', 'Night period-product gap (pp)', 2, 'day_product_gap', 'Day period-product gap (pp)', 2,
  'day_traffic_share', 'Traffic share of daytime list (\\%)', 1, 'night_traffic_share', 'Traffic share of night list (\\%)', 1)
tex(data.frame(Quantity = q$label, All = map2_chr(q$col, q$k, ~fmt(fs[[.x]], .y)), Entrance = map2_chr(q$col, q$k, ~rng(min(oe[[.x]]), max(oe[[.x]]), .y)),
  Date = map2_chr(q$col, q$k, ~rng(min(od[[.x]]), max(od[[.x]]), .y))) %>% setNames(c('Quantity', 'All entrances', 'One entrance omitted', 'One date omitted')), file.path(tab, 'si/omissions.tex'))
write_csv(hr, file.path(out, 'hourly_summary.csv')); write_csv(filter(rules, budget == 10), file.path(out, 'selection_ten.csv')); write_csv(de, file.path(out, 'decomposition.csv'))
# Reference tables copied from previously verified original outputs, never retyped.
old <- 'intermediate/support_tables'
for (n in c('si_spatial.tex', 'si_beyond_loudness.tex')) stopifnot(file.copy(file.path(old, n), file.path(tab, 'si', n), overwrite = TRUE))
# Validation at the measured unit; no citywide hourly-validation claim.
vpath <- file.path(e, 'analysis/33_sound_event_validation__summary.csv'); cvpath <- file.path(e, 'analysis/21_domain_cv__model_performance.csv')
v <- read_csv(vpath, show_col_types = FALSE) %>% filter(subset %in% c('all', 'night')) %>% arrange(level, task, subset)
tex(data.frame(Resolution = ifelse(v$level == 'group', 'Four groups', 'Twelve categories'), Label = ifelse(v$task == 'pri', 'Primary', 'Secondary'), Subset = ifelse(v$subset == 'all', 'All clips', 'Night'), N = v$N_valid_majority,
  Accuracy = sprintf('%.3f', v$accuracy), Kappa = sprintf('%.3f', v$kappa)), file.path(tab, 'si/validation.tex'), 'lllrrr')
cv <- read_csv(cvpath, show_col_types = FALSE)
mod <- c(intercept_only = 'Intercept', environment = 'Context', source = 'Source', environment_source = 'Context + source')
periodlab <- c(laeq_24h_leq_db_a = '24 h', lday_leq_db_a = 'Day', lnight_leq_db_a = 'Night', lpeak_leq_db_a = 'Peak', linter_peak_leq_db_a = 'Inter-peak')
cv <- cv %>% filter(model %in% names(mod), outcome %in% names(periodlab)) %>% arrange(match(outcome, names(periodlab)), match(model, names(mod)))
tex(data.frame(Period = unname(periodlab[cv$outcome]), Model = unname(mod[cv$model]), Predictors = cv$p, Adjusted = sprintf('%.3f', ifelse(abs(cv$adj_r2) < 1e-12, 0, cv$adj_r2)), RMSE = paste0(sprintf('%.2f', cv$cv_rmse_mean), ' $\\pm$ ', sprintf('%.2f', cv$cv_rmse_sd))), file.path(tab, 'si/context.tex'), 'llrrr')
# Every display has a file-level provenance record for the generated displays.
inputs <- c(file.path(e, 'data/clean', d$input_hashes$file), roadpath, used, cvpath, vpath, file.path(old, c('si_spatial.tex', 'si_beyond_loudness.tex')), script)
write_csv(d$input_hashes, file.path(out, 'clean_data_sha256.csv'))
write_csv(tibble(path = sub(paste0(project, '/'), '', inputs, fixed = TRUE), sha256 = vapply(inputs, function(f) digest::digest(file = f, algo = 'sha256'), character(1))), file.path(out, 'display_inputs_sha256.csv'))
writeLines(capture.output(sessionInfo()), file.path(out, 'session_info.txt'))
cat('Built seven main figures, two SI figures and table fragments.\n')
