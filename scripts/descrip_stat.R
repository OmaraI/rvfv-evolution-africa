
## Descriptive statistics for all countries

# Load libraries
library(dplyr)
library(stats)
library(rstatix)   # for pairwise tests + multiple testing correction
library(readr)

# Get WD
getwd()

# Set the WD
setwd("~/Desktop/Manuscripts/Data/R-Work/")

## Read in your dataset
data <- read_csv("Combined-Mutation.csv")

## Inspect structure
head(data)


## Descriptive statistics by country

summary_stats <- data %>%
  group_by(Country) %>%
  summarise(
    n = n(),
    mean = mean(`Total Protein Mutations`, na.rm = TRUE),
    sd = sd(`Total Protein Mutations`, na.rm = TRUE),
    median = median(`Total Protein Mutations`, na.rm = TRUE),
    q1 = quantile(`Total Protein Mutations`, 0.25, na.rm = TRUE),
    q3 = quantile(`Total Protein Mutations`, 0.75, na.rm = TRUE)
  ) %>%
  arrange(desc(median))

print(summary_stats)


## Kruskal-Wallis test
kw_test <- kruskal_test(data, `Total Protein Mutations` ~ Country)
print(kw_test)


## Pairwise Wilcoxon tests

pairwise_wilcox <- data %>%
  pairwise_wilcox_test(
    `Total Protein Mutations` ~ Country,
    p.adjust.method = "BH"   # Benjamini-Hochberg correction
  )

print(pairwise_wilcox)


## Optional: Export results
write_csv(summary_stats, "RVFV_mutation_summary_by_country.csv")
write_csv(pairwise_wilcox, "RVFV_pairwise_wilcoxon_results.csv")
