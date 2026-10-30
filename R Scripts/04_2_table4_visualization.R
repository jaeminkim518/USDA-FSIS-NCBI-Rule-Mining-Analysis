## 04_2_table4_plots.R
## Plots of the Table 4 files (extension of 04_table4.R)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/USDA/Rule Mining Work")


# Read the 5 Table 4 files and put them together
beef <- read_csv("Output/table4_Beef.csv")
chicken <- read_csv("Output/table4_Chicken.csv")
pork <- read_csv("Output/table4_Pork.csv")
turkey <- read_csv("Output/table4_Turkey.csv")
clinical <- read_csv("Output/table4_Clinical.csv")

table4 <- bind_rows(beef, chicken, pork, turkey, clinical)

# Remove the Mean and SD rows, then turn Year back into a number
table4 <- table4 %>% filter(Year != "Mean" & Year != "SD")
table4$Year <- as.numeric(table4$Year)


# One row per class, so ggplot can draw one line per class
# (macrolide and colistin left out, they're near 0% everywhere)
long <- table4 %>%
  select(specimen_type, sample_type, Year, AMINOGLYCOSIDE, `BETA-LACTAM`,
         `FOLATE-PATHWAY-INHIBITOR`, PHENICOL, QUINOLONE, TETRACYCLINE) %>%
  pivot_longer(cols = AMINOGLYCOSIDE:TETRACYCLINE, names_to = "class", values_to = "percent")

head(long)


# Resistance by class: one plot per sample type, one panel per specimen type
for (st in c("Beef", "Chicken", "Pork", "Turkey", "Clinical")) {
  st_data <- long %>% filter(sample_type == st)
  
  p <- ggplot(st_data, aes(x = Year, y = percent, color = class)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ specimen_type) +
    scale_x_continuous(breaks = seq(2015, 2026, 2)) +   # whole years on the x-axis
    labs(y = "% of isolates", title = paste("Resistance by class -", st))
  
  ggsave(paste0("Output/table4_plot_", st, ".png"), p, width = 10, height = 4)
}


# MDR by year, all specimen types on one plot
ggplot(table4, aes(x = Year, y = MDR, color = specimen_type)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ sample_type) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(y = "% MDR (3+ classes)", title = "MDR by year")

ggsave("Output/mdr_by_year.png", width = 10, height = 6)

