## Figure 2: Distribution of protein-coding mutations per isolate across countries

# Set the working directory
setwd("~/Desktop/Manuscripts/Data/R-Work/") # Adjust accordingly

# Load the libraries
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# Read CSVs
L_NH <- read_csv("L-NH-Africa.csv") %>% mutate(Host = "Non-human")
M_NH <- read_csv("M-NH-Africa.csv") %>% mutate(Host = "Non-human")
S_H  <- read_csv("S_H.csv")         %>% mutate(Host = "Human")
S_NH <- read_csv("S_NH.csv")        %>% mutate(Host = "Non-human")
L_H  <- read_csv("L-H-Africa.csv")  %>% mutate(Host = "Human")
M_H  <- read_csv("M-H-Africa.csv")  %>% mutate(Host = "Human")

# Combine & clean (assumes L_H, L_NH, M_H, M_NH, S_H, S_NH)
mutation_data <- bind_rows(L_H, L_NH, M_H, M_NH, S_H, S_NH) %>%
  mutate(
    Country = str_trim(Country),
    Host    = str_trim(Host),
    Segment = str_to_upper(str_trim(Segment)),
    `Total Protein Mutations` = suppressWarnings(as.numeric(`Total Protein Mutations`))
  ) %>%
  filter(!is.na(`Total Protein Mutations`), !is.na(Country))

# Save/Load
write_csv(mutation_data, "Combined-Mutation.csv")

## Box-plot:
# Boxplot by Country (y-axis starts at 0, countries alphabetical)

p_box <- ggplot(
  mutation_data,
  aes(
    x = factor(Country, levels = sort(unique(Country))),  # <- alphabetical order
    y = `Total Protein Mutations`,
    fill = Country
  )
) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5, alpha = 0.8) +
  labs(
    title = "Protein Coding Mutation Burden per Isolate, by Country",
    x = "Country",
    y = "Total Protein Mutations"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 12)),
    axis.title  = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain"),  # alphabet labels
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.line = element_line(color = "black", linewidth = 0.3)
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))

print(p_box)
ggsave("Boxplot_Mutation_By_Country.pdf", p_box, width = 10, height = 6, dpi = 300)
