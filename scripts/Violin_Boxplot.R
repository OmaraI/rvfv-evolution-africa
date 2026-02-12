## Violin + Boxplot: Distribution of Protein-Coding Mutations by Segment and Host 

# Libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(janitor)
})

# Load combined dataset
mutation_data <- read_csv("Combined-Mutation.csv", show_col_types = FALSE) %>%
  clean_names()

# Expect (after clean_names): country, host, segment, total_protein_mutations, etc.

# Clean host/segment + numeric mutations
mutation_data <- mutation_data %>%
  mutate(
    segment = str_to_upper(str_squish(segment)),
    host_raw = str_to_lower(str_squish(host)),
    host = case_when(
      host_raw %in% c("human", "homo sapiens", "homo_sapiens") ~ "Human",
      host_raw %in% c("non-human", "non human", "nonhuman", "animal", "livestock", "wildlife", "vector") ~ "Non-human",
      TRUE ~ "Unknown"
    ),
    total_protein_mutations = suppressWarnings(as.numeric(total_protein_mutations))
  ) %>%
  filter(
    segment %in% c("L", "M", "S"),
    host %in% c("Human", "Non-human"),
    !is.na(total_protein_mutations)
  )

# Summary table
summary_tbl <- mutation_data %>%
  group_by(segment, host) %>%
  summarise(
    n      = n(),
    mean   = mean(total_protein_mutations, na.rm = TRUE),
    median = median(total_protein_mutations, na.rm = TRUE),
    IQR    = IQR(total_protein_mutations, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tbl)

# Plot: Violin + Box + Jitter
plot_vb <- ggplot(
  mutation_data,
  aes(x = segment, y = total_protein_mutations, fill = host)
) +
  geom_violin(
    position = position_dodge(width = 0.8),
    trim = FALSE,
    alpha = 0.55
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.9,
    linewidth = 0.3
  ) +
  geom_jitter(
    position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.15),
    size = 0.8,
    alpha = 0.35,
    show.legend = FALSE,
    color = "grey35"
  ) +
  labs(
    title = "Distribution of Protein-Coding Mutations by Segment and Host",
    x = "Genome Segment",
    y = "Total Protein Mutations"
  ) +
  # Keep your palette exactly (Human blue, Non-human orange)
  scale_fill_manual(
    values = c("Human" = "#1f77b4", "Non-human" = "#ff7f0e"),
    name = "Host"
  ) +
  # Ensure segment order L-M-S
  scale_x_discrete(limits = c("L", "M", "S")) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 10)),
    axis.title   = element_text(face = "bold", size = 12),
    axis.text    = element_text(size = 11),
    legend.title = element_text(face = "bold", size = 11),
    legend.text  = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.line    = element_line(color = "black", linewidth = 0.3)
  )

print(plot_vb)

# Save
ggsave("Violin-Boxplot_segment_host.pdf", plot = plot_vb,
       width = 9, height = 6, units = "in", dpi = 300)

ggsave("Violin-Boxplot_segment_host.png", plot = plot_vb,
       width = 9, height = 6, units = "in", dpi = 600, bg = "white")
