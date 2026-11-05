## 07_rule_mining.R
## Class-level association rule mining for chicken

library(readr)
library(dplyr)
library(arules)

setwd("~/USDA/Rule Mining Work")

isolate_classes <- read_csv("Output/isolate_classes.csv")
isolates <- read_csv("Output/isolates_clean.csv",
                     col_select = c(BioSample, specimen_type, sample_type, Year))

narms_classes <- c("AMINOGLYCOSIDE", "BETA-LACTAM", "FOLATE-PATHWAY-INHIBITOR", "MACROLIDE",
                   "PHENICOL", "COLISTIN", "QUINOLONE", "TETRACYCLINE")

dat <- isolate_classes %>%
  left_join(isolates, by = "BioSample") %>%
  filter(sample_type == "Chicken",
         specimen_type %in% c("cecal", "non-cecal", "retail"))

dat$group <- paste(dat$specimen_type, dat$Year)

mine_group <- function(df) {
  m <- as.matrix(df[, narms_classes])
  mode(m) <- "logical"
  trans <- as(m, "transactions")
  
  n <- nrow(df)
  rules <- apriori(trans,
                   parameter = list(support = 1/n, confidence = 0,
                                    target = "rules", minlen = 2),
                   control = list(verbose = FALSE))
  
  if (length(rules) == 0) return(rules)
  
  q <- interestMeasure(rules, c("cosine", "jaccard", "kulczynski"))
  quality(rules) <- cbind(quality(rules), q)
  subset(rules, cosine > 0.5 & kulczynski > 0.5 & jaccard > 0 & support > 0)
}

# Mine every group with at least 20 isolates
all_rules <- list()
for (g in unique(dat$group)) {
  gdf <- dat %>% filter(group == g)
  if (nrow(gdf) < 20) next
  r <- mine_group(gdf)
  if (length(r) == 0) next
  rdf <- DATAFRAME(r)
  rdf$group <- g
  rdf$n_isolates <- nrow(gdf)
  all_rules[[g]] <- rdf
}

rules_out <- bind_rows(all_rules)
write_csv(rules_out, "Output/rules_chicken_class.csv")