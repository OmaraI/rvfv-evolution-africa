## Spatiotemporal distribution of RVFV protein-coding mutations in Africa

# Packages required
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)   # for side-by-side layout

# 1) Set the working directory
setwd("~/Desktop/Manuscripts/Data/R-Work/")

# Load data
mutation_data <- read_csv("Combined-Mutation.csv")

# Ensure numeric types and basic filtering
mutation_data <- mutation_data %>%
  mutate(
    Year = suppressWarnings(as.integer(`Year of Sample Collection`)),
    `Total Protein Mutations` = suppressWarnings(as.numeric(`Total Protein Mutations`))
  ) %>%
  filter(!is.na(Country), !is.na(Year), !is.na(`Total Protein Mutations`))


# HEATMAP (left)
# Geographic & temporal distribution of mean mutation burden
heatmap_data <- mutation_data %>%
  group_by(Country, Year) %>%
  summarise(Mean_Mutations = mean(`Total Protein Mutations`, na.rm = TRUE), .groups = "drop")

# Order countries from LOW to HIGH average burden (across all years)
country_order <- heatmap_data %>%
  group_by(Country) %>%
  summarise(Total = mean(Mean_Mutations, na.rm = TRUE), .groups = "drop") %>%
  arrange(Total) %>%
  pull(Country)

# Enforce the same country order for ALL downstream plots
mutation_data <- mutation_data %>%
  mutate(Country = factor(Country, levels = country_order))

heatmap_data$Country <- factor(heatmap_data$Country, levels = country_order)

heatmap_plot <- ggplot(heatmap_data, aes(x = Year, y = Country, fill = Mean_Mutations)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "#440154", high = "#FDE725", name = "Avg Mutations") +
  labs(
    title = "1a",
    #subtitle = "Mean total protein mutations per country per year",
    x = "Year of Sample Collection", y = "Country"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.0, face = "bold", size = 14, margin = margin(b = 4)),
    plot.subtitle= element_text(hjust = 0.0, size = 11, margin = margin(b = 10)),
    
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    
    axis.text.x  = element_text(size = 16, angle = 40, hjust = 1),
    axis.text.y  = element_text(
      size = 16,                # Larger country labels
      margin = margin(r = 8),   # Pushes countries away from heatmap
      lineheight = 1.15         # Adds vertical spacing for tall labels
    ),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.2),
    
    plot.margin  = margin(12, 12, 12, 12)
  )

# -------------------------------------------------
# WATERFALL (right)
#    Country-level YoY changes and start/end totals
# -------------------------------------------------
mutation_trends <- mutation_data %>%
  group_by(Country, Year) %>%
  summarise(Mean_Mutations = mean(`Total Protein Mutations`, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(Year)) %>%
  arrange(Country, Year)

build_waterfall <- function(df_country) {
  df_country <- df_country %>% arrange(Year)
  if (nrow(df_country) == 0) return(tibble())
  
  start_total <- df_country$Mean_Mutations[1]
  end_total   <- df_country$Mean_Mutations[nrow(df_country)]
  
  change_rows <- tibble(
    label = as.character(df_country$Year[-1]),
    value = df_country$Mean_Mutations[-1] - df_country$Mean_Mutations[-nrow(df_country)],
    type  = "change"
  )
  
  wf <- bind_rows(
    tibble(label = paste0(df_country$Year[1], " (Start)"),
           value = start_total, type = "total_start"),
    change_rows,
    tibble(label = paste0(df_country$Year[nrow(df_country)], " (End)"),
           value = end_total, type = "total_end")
  ) %>%
    mutate(x = row_number(),
           running = ifelse(type == "total_start", value, NA_real_))
  
  # running total
  for (i in seq_len(nrow(wf))) {
    if (i == 1) next
    if (wf$type[i] == "change") {
      wf$running[i] <- wf$running[i-1] + wf$value[i]
    } else if (wf$type[i] == "total_end") {
      wf$running[i] <- wf$value[i]
    }
  }
  
  wf %>%
    mutate(
      prev_running = lag(running, default = 0),
      ymin = case_when(
        type %in% c("total_start", "total_end") ~ 0,
        value >= 0 ~ prev_running,
        TRUE ~ running
      ),
      ymax = case_when(
        type %in% c("total_start", "total_end") ~ value,
        value >= 0 ~ running,
        TRUE ~ prev_running
      ),
      fill_type = case_when(
        type == "change" & value >= 0 ~ "Increase",
        type == "change" & value < 0  ~ "Decrease",
        type %in% c("total_start", "total_end") ~ "Total"
      )
    )
}

waterfall_df <- mutation_trends %>%
  group_by(Country) %>%
  group_modify(~ build_waterfall(.x)) %>%
  ungroup()

# Reorder countries for waterfall ONLY: highest -> lowest
waterfall_order <- heatmap_data %>%
  group_by(Country) %>%
  summarise(Total = mean(Mean_Mutations, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%     # HIGH -> LOW
  pull(Country)

waterfall_df <- waterfall_df %>%
  mutate(Country = factor(Country, levels = waterfall_order))


y_max <- ceiling(max(waterfall_df$ymax, na.rm = TRUE) / 5) * 5

waterfall_plot <- ggplot(waterfall_df) +
  geom_rect(aes(xmin = x - 0.45, xmax = x + 0.45, ymin = ymin, ymax = ymax, fill = fill_type),
            color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(aes(x = x, y = pmax(ymin, ymax),
                label = ifelse(type == "change", sprintf("%+.1f", value), sprintf("%.1f", ymax))),
            vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = c("Increase" = "#0072B2", "Decrease" = "#E69F00", "Total" = "#6D6D6D")) +
  scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, by = 5)) +
  labs(
    title = "1b",
    #subtitle = "Start = first year total; steps = YoY change; End = last year total",
    x = NULL, y = "Mean Total Protein Mutations", fill = "Change Type"
  ) +
  facet_wrap(~ Country, scales = "free_x") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.0, face = "bold", size = 14, margin = margin(b = 2)),
    plot.subtitle= element_text(hjust = 0.0, size = 10, margin = margin(b = 8)),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line    = element_line(color = "black", linewidth = 0.4),
    legend.title = element_text(face = "bold"),
    strip.text   = element_text(face = "bold", size = 12),
    plot.margin  = margin(10, 10, 10, 10)
  )

# ------------------------------------------
# COMBINE (heatmap left, waterfall right)
# ------------------------------------------
# widths control the space each plot takes; tweak as you like (e.g., c(1, 1.25))
combined <- heatmap_plot + waterfall_plot +
  plot_layout(widths = c(2, 2)) &
  theme(legend.position = "right")  # keep legends right-aligned

print(combined)

# -----------------------------------------
# Save outputs
# -----------------------------------------
ggsave("Combined_Heatmap_Waterfall.pdf", plot = combined, width = 20, height = 8.5, dpi = 300)
