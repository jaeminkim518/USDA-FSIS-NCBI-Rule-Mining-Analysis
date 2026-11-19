## 08_rule_comparison.R
## Comparing chicken rulesets across years & across stages
## Calculate roportion captured & cumulative proportion captured (CPC) following Glass et al.

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("~/USDA/Rule Mining Work")

rules <- read_csv("Output/rules_chicken_class_serotype_test.csv")

rules$rule <- paste(rules$LHS, "=>", rules$RHS)
rules$stage <- word(rules$group, 1)
rules$Year <- as.numeric(word(rules$group, 2))

# (Glass et al. Eq. 5)
prop_captured <- function(r1, r2) {
  length(intersect(r1, r2)) / length(r1)
}

# rules in both / rules in either (project plan)
overlap <- function(r1, r2) {
  length(intersect(r1, r2)) / length(union(r1, r2))
}


# Across years (Glass et al. Eq. 6)
time_comp <- data.frame()

for (s in unique(rules$subset)) {
  for (st in unique(rules$stage)) {
    d <- rules %>% filter(subset == s, stage == st)
    years <- sort(unique(d$Year))
    if (length(years) < 2) next
    
    for (i in 1:(length(years) - 1)) {
      r1 <- d$rule[d$Year == years[i]]
      r2 <- d$rule[d$Year == years[i + 1]]
      pc <- prop_captured(r1, r2)
      
      if (i == 1) {
        cpc <- pc
      } else {
        cpc <- (1 / i) * ((i - 1) * cpc + pc)
      }
      
      time_comp <- rbind(time_comp,
                         data.frame(subset = s, stage = st,
                                    year1 = years[i], year2 = years[i + 1],
                                    n_rules1 = length(r1), n_rules2 = length(r2),
                                    prop_captured = pc, cpc = cpc,
                                    overlap = overlap(r1, r2)))
    }
  }
}

time_comp
write_csv(time_comp, "Output/rule_comparison_time.csv")


# Across stages (Glass et al. Eq. 7)
stage_pairs <- list(c("cecal", "non-cecal"), c("non-cecal", "retail"), c("cecal", "retail"))

stage_comp <- data.frame()

for (s in unique(rules$subset)) {
  for (p in stage_pairs) {
    d1 <- rules %>% filter(subset == s, stage == p[1])
    d2 <- rules %>% filter(subset == s, stage == p[2])
    years <- sort(intersect(d1$Year, d2$Year))
    if (length(years) == 0) next
    
    for (i in 1:length(years)) {
      r1 <- d1$rule[d1$Year == years[i]]
      r2 <- d2$rule[d2$Year == years[i]]
      pc_12 <- prop_captured(r1, r2)
      pc_21 <- prop_captured(r2, r1)
      
      if (i == 1) {
        cpc_12 <- pc_12
        cpc_21 <- pc_21
      } else {
        cpc_12 <- (1 / i) * ((i - 1) * cpc_12 + pc_12)
        cpc_21 <- (1 / i) * ((i - 1) * cpc_21 + pc_21)
      }
      
      stage_comp <- rbind(stage_comp,
                          data.frame(subset = s, Year = years[i],
                                     pair = paste(p[1], "vs", p[2]),
                                     direction = c(paste(p[1], "in", p[2]),
                                                   paste(p[2], "in", p[1])),
                                     type = c("Upstream rules found downstream",
                                              "Downstream rules found upstream"),
                                     prop_captured = c(pc_12, pc_21),
                                     cpc = c(cpc_12, cpc_21),
                                     overlap = overlap(r1, r2)))
    }
  }
}

stage_comp
write_csv(stage_comp, "Output/rule_comparison_stage.csv")


# Plots
stage_colors <- c("cecal" = "#C8A165", "non-cecal" = "#B4645C", "retail" = "#5B7FA6")

# year to year
ggplot(time_comp, aes(x = year2, y = prop_captured, color = stage)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ subset) +
  scale_color_manual(values = stage_colors) +
  scale_x_continuous(breaks = seq(2016, 2026, 2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Year (compared with previous year)",
       y = "Share of previous year's rules still found",
       color = "Stage", title = "Rule stability year to year (chicken)")

ggsave("Output/Figures/rule_comparison_time.png", width = 10, height = 6)

# cumulative (Glass et al. Fig. 2)
ggplot(time_comp, aes(x = year2, y = cpc, color = stage)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ subset) +
  scale_color_manual(values = stage_colors) +
  scale_x_continuous(breaks = seq(2016, 2026, 2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Year", y = "Cumulative proportion captured",
       color = "Stage", title = "Cumulative rule stability over time (chicken)")

ggsave("Output/Figures/rule_comparison_time_cpc.png", width = 10, height = 6)

# between stages (Glass et al. Fig. 3)
stage_comp$pair <- factor(stage_comp$pair,
                          levels = c("cecal vs non-cecal", "non-cecal vs retail", "cecal vs retail"))

ggplot(stage_comp, aes(x = Year, y = cpc, color = type)) +
  geom_line() +
  geom_point(size = 0.8) +
  facet_grid(pair ~ subset) +
  scale_color_manual(values = c("Upstream rules found downstream" = "#B4645C",
                                "Downstream rules found upstream" = "#5B7FA6")) +
  scale_x_continuous(breaks = seq(2015, 2026, 3)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Cumulative proportion captured", color = NULL,
       title = "Rules shared between chicken stages",
       caption = "Order: cecal -> non-cecal -> retail") +
  theme(legend.position = "bottom")

ggsave("Output/Figures/rule_comparison_stage.png", width = 12, height = 8)