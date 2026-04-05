
### Mutations assoiciated with viral evolution & host adaptation — Figure 6

# ------------------------------
# 1) Libraries & Data Setup
# ------------------------------
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(scales)
library(ggrepel)


# ------------------------------
# 1) Data Prep 
# ------------------------------
metadata <- read_csv("~/Desktop/Manuscripts/Data/R-Work/Combined-Mutation.csv", show_col_types = FALSE)

metadata_long <- metadata %>%
  mutate(Year = suppressWarnings(as.integer(`Year of Sample Collection`)),
         Mutation = str_split(`Protein Mutations`, ";|,")) %>%
  unnest(Mutation) %>%
  mutate(Mutation = str_trim(Mutation)) %>%
  filter(!is.na(Mutation), Mutation != "")

mutations <- c("N277S","N277D","S278N","I442S","I442V","V659A","N133S")
mutation_hits <- metadata_long %>% filter(Mutation %in% mutations)

mutation_summary <- mutation_hits %>%
  group_by(Mutation, Country, Year) %>%
  summarise(n_sequences = n(), .groups = "drop") %>%
  filter(!is.na(Year))

# Factors and Colors
country_order <- mutation_summary %>% group_by(Country) %>% summarise(t = sum(n_sequences)) %>% arrange(desc(t)) %>% pull(Country)
mutation_summary <- mutation_summary %>% mutate(Country = factor(Country, levels = country_order), Mutation = factor(Mutation, levels = mutations))
cb_palette <- c("N277S"="#1b9e77","N277D"="#d95f02","S278N"="#7570b3","I442S"="#e7298a","I442V"="#66a61e","V659A"="#e6ab02","N133S"="#a6761d")

# ======================
#### PANEL A — AFRICA MAP 
# ======================
africa <- ne_countries(scale = "medium", returnclass = "sf") %>% filter(region_un == "Africa")
africa_data <- africa %>% left_join(mutation_summary %>% group_by(Country) %>% summarise(total = sum(n_sequences)), by = c("name" = "Country"))

## Filter for active countries
active_countries <- africa_data %>%
  filter(total > 0)

map_plot <- ggplot(africa_data) +
  geom_sf(aes(fill = total), color = "grey30", linewidth = 0.1) +
  scale_fill_gradientn(
    colors = c("#4d0000", "#b30000", "#ff6600", "#ffcc00", "#ffff99"),
    trans  = "pseudo_log",
    breaks = c(2, 5, 10, 20, 50, 100),
    labels = label_number(accuracy = 1),
    na.value = "#f0f0f0" 
  ) +
  geom_text_repel(
    data = active_countries,
    aes(label = name, geometry = geometry),
    stat = "sf_coordinates",
    size = 5,             
    fontface = "bold",
    color = "black",
    bg.color = "white",
    bg.r = 0.1,
    box.padding = 0.3,
    max.overlaps = Inf
  ) +
  coord_sf(expand = FALSE) + 
  labs(title = "(a)", fill = "Total\nSequences") +
  theme_void() +
  theme(
    plot.title = element_text(size = 32, face = "bold"),
    legend.position = c(0.12, 0.25), # Tucked into bottom left
    legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    # Aggressive margin to let the map use every pixel
    plot.margin = margin(r = -110, l = 10, t = 10, b = 10) 
  ) +
  guides(fill = guide_colorbar(barwidth = 1.8, barheight = 8, frame.colour = "black"))

# =================------
## PANEL B — BAR PLOT 
# =================------
bar_plot <- ggplot(mutation_summary, 
                   aes(x = Country, y = n_sequences, fill = Mutation)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = cb_palette) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1), 
    breaks = c(1, 10, 100),
    labels = label_number()
  ) +
  labs(title = "(b)", y = "Sequences", x = NULL) +
  theme_minimal(base_size = 16) + 
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold", color = "black"),
    plot.margin = margin(l = -110, r = 10, t = 10, b = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# =================------------
# PANEL C — TEMPORAL 
# =================------------
timeline_df <- mutation_summary %>% group_by(Year, Mutation) %>% summarise(n = sum(n_sequences), .groups = "drop")

timeline_plot <- ggplot(timeline_df, aes(x = Year, y = n, color = Mutation)) +
  geom_line(linewidth = 2) + 
  geom_point(size = 4.5) +
  scale_color_manual(values = cb_palette) +
  scale_x_continuous(limits = c(1940, 2020), breaks = seq(1940, 2020, 10)) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1), 
    breaks = c(1, 10, 100),
    labels = label_number()
  ) +
  labs(title = "(c)", x = "Year", y = "Sequences") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold", color = "black"),
    plot.margin = margin(l = -110, r = 10, t = 10, b = 10)
  )

# =================--------------
# ASSEMBLE FIGURE 
# =================--------------
# 'AAAAAA' (6 units) vs 'BB' (2 units) makes the map much larger
design <- "
  AAAAAABBB
  AAAAAACCC
"

fig_final <- wrap_plots(A = map_plot, B = bar_plot, C = timeline_plot, design = design)

# -------------
# Final Save
# -------------
# Increased width to 28 to ensure the legend on the far right isn't squeezed
ggsave("Figure6.png", plot = fig_final, width = 28, height = 14, dpi = 300)

print(fig_final)
