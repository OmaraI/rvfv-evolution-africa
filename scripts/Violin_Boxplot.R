## Violin-Boxplot: Distribution of Protein-Coding Mutations by host and segment ###

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(janitor)
  library(broom)
})

setwd("~/Desktop/Manuscripts/Data/R-Work/") # Update accordingly

# Read CSVs
L_H  <- read_csv("L-H-Africa.csv") %>% mutate(Host = "human")
L_NH <- read_csv("L-NH-Africa.csv")
M_H  <- read_csv("M-H-Africa.csv") %>% mutate(Host = "human")
M_NH <- read_csv("M-NH-Africa.csv")
S_H  <- read_csv("S_H.csv")  %>% mutate(Host = "human")
S_NH <- read_csv("S_NH.csv") %>% mutate(Host = "Non-human")

# Combine all datasets
mutation_data <- bind_rows(L_H, L_NH, M_H, M_NH, S_H, S_NH)

# Save the mutation data
write_csv(mutation_data, "Combined-Mutation.csv")

# Load the combined mutation dataset
mutation_data <- read_csv("Combined-Mutation.csv") 
# Load combined dataset
mutation_data <- read_csv("Combined-Mutation.csv", show_col_types = FALSE) %>%
  clean_names()

# Clean host/segment + numeric mutations
mutation_data <- mutation_data %>%
  mutate(
    segment = str_to_upper(str_squish(segment)),
    host_raw = str_to_lower(str_squish(host)),
    host = case_when(
      host_raw %in% c("human", "homo sapiens", "homo_sapiens") ~ "Human",
      host_raw %in% c("non-human", "non human", "nonhuman", "animal",
                      "livestock", "wildlife", "vector") ~ "Non-human",
      TRUE ~ "Unknown"
    ),
    total_protein_mutations = suppressWarnings(as.numeric(total_protein_mutations))
  ) %>%
  filter(
    segment %in% c("L", "M", "S"),
    host %in% c("Human", "Non-human"),
    !is.na(total_protein_mutations),
    total_protein_mutations >= 0
  )

# Summary table
summary_tbl <- mutation_data %>%
  group_by(segment, host) %>%
  summarise(
    n = n(),
    mean = mean(total_protein_mutations, na.rm = TRUE),
    median = median(total_protein_mutations, na.rm = TRUE),
    IQR = IQR(total_protein_mutations, na.rm = TRUE),
    q1 = quantile(total_protein_mutations, 0.25, na.rm = TRUE),
    q3 = quantile(total_protein_mutations, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tbl)
write_csv(summary_tbl, "Figure3_segment_host_summary_stats.csv")

# Wilcoxon tests: Human vs Non-human within each segment
wilcox_tbl <- mutation_data %>%
  group_by(segment) %>%
  group_modify(~{
    test <- wilcox.test(
      total_protein_mutations ~ host,
      data = .x,
      conf.int = TRUE,
      exact = FALSE
    )
    tibble(
      p_value = test$p.value,
      median_difference = unname(test$estimate),
      conf_low = test$conf.int[1],
      conf_high = test$conf.int[2]
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_value_formatted = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

print(wilcox_tbl)
write_csv(wilcox_tbl, "Figure3_wilcoxon_host_comparisons.csv")

# Plot: lighter violin + clearer boxplot + jitter
plot_vb <- ggplot(
  mutation_data,
  aes(x = segment, y = total_protein_mutations, fill = host)
) +
  geom_violin(
    position = position_dodge(width = 0.8),
    trim = TRUE,
    alpha = 0.25,
    color = "grey35",
    linewidth = 0.25
  ) +
  geom_boxplot(
    aes(color = host),
    position = position_dodge(width = 0.8),
    width = 0.16,
    outlier.shape = NA,
    alpha = 1,
    fill = "white",
    linewidth = 0.65
  ) +
  geom_jitter(
    position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.12),
    size = 0.8,
    alpha = 0.22,
    show.legend = FALSE,
    color = "grey30"
  ) +
  labs(
    #title = "Distribution of Protein-Coding Mutations by Segment and Host",
    x = "Genome Segment",
    y = "Total Protein Mutations"
  ) +
  scale_fill_manual(
    values = c("Human" = "#9ecae1", "Non-human" = "#fdd0a2"),
    name = "Host"
  ) +
  scale_color_manual(
    values = c("Human" = "#1f77b4", "Non-human" = "#ff7f0e"),
    guide = "none"
  ) +
  scale_x_discrete(limits = c("L", "M", "S")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.line = element_line(color = "black", linewidth = 0.3)
  )

print(plot_vb)

ggsave("Violin-Boxplot_segment_host_revised.pdf", plot = plot_vb,
       width = 9, height = 6, units = "in", dpi = 300)

ggsave("Violin-Boxplot_segment_host_revised.png", plot = plot_vb,
       width = 9, height = 6, units = "in", dpi = 600, bg = "white")
