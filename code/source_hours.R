suppressPackageStartupMessages(library(tidyverse))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- "intermediate"
source("code/shared/helpers.R")
source("code/shared/plot_style.R")
d <- load_p11(root)
out <- file.path(root, "analysis_38_priority_sources/outcome")
s <- p11_scores(d$panel) %>% group_by(period) %>% mutate(rank = rank(-raw50_joint, ties.method = "first")) %>% ungroup()
ids <- s %>% filter(period %in% c("day", "night"), rank <= 10) %>% select(station_key, period, rank) %>%
  pivot_wider(names_from = period, values_from = rank) %>% arrange(is.na(day), day, night) %>%
  mutate(label = paste0(station_key, "  [", ifelse(is.na(day), "–", day), "/", ifelse(is.na(night), "–", night), "]"))
p <- d$panel %>% inner_join(ids, by = "station_key") %>% mutate(label = factor(label, levels = rev(ids$label)),
  source_group = factor(source_group, levels = c("traffic", "human", "natural", "mechanical")))
g <- ggplot(p, aes(hour, label, fill = source_group)) + geom_tile(colour = "white", linewidth = .1) +
  geom_vline(xintercept = c(6.5, 22.5), linewidth = .4) +
  scale_fill_manual(values = c(traffic = "#0072B2", human = "#E69F00", natural = "#009E73", mechanical = "#CC79A7"),
    labels = c("Traffic", "Human", "Natural", "Mechanical")) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 23), expand = expansion(add = .5)) +
  labs(x = "Clock hour", y = "Entrance [day rank / night rank]") + p11_theme() +
  theme(axis.ticks.y = element_blank(), legend.text = element_text(size = 9))
p11_png(g, file.path(out, "priority_source_hours.png"))
