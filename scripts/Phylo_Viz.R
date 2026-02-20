# Phylogenetic Visualization and Metadata Annotation [RVFV M segment] 
[Tree + Region + Lineage] + bootstraps

	•	Mapping bootstrap values
	•	Coloring tips by Host
	•	Adding Region and Lineage heatmaps

## Load required libraries
library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(phytools)

## Set working directory (adjust as needed)
setwd("~/IsaacOmara/Desktop/Manuscripts/Data/R-Work/Nexus-Files/M/Bootstrap/") # Adjust

# Paths (adjust if needed)
iqtree_file <- "M_all.treefile"
timetree_file <- "../treetime_M/timetree.nexus"

# Read trees
t_iq <- read.tree(iqtree_file)            # has bootstrap supports in node.label
t_tt <- read.nexus(timetree_file)         # time-scaled (likely no supports)

# Sanity: same tips?
stopifnot(setequal(t_iq$tip.label, t_tt$tip.label))

# Optional: root both the same way if they differ (often not needed)
# t_iq <- midpoint.root(t_iq)
# t_tt <- midpoint.root(t_tt)

# Parse support from IQ-TREE node labels (handles "95" or "95/0.98" formats)
parse_support <- function(x) {
  if (is.null(x)) return(rep(NA_real_, 0))
  y <- sub("^\\s*([0-9]+).*$", "\\1", x)       # take leading digits
  suppressWarnings(as.numeric(y))
}
iq_bs <- parse_support(t_iq$node.label)

# Match internal nodes between timetree (target) and iqtree (source)
map <- matchNodes(t_tt, t_iq)  # columns: t_tt node, t_iq node (both as node numbers)

# Prepare vector for timetree supports
tt_bs <- rep(NA_real_, t_tt$Nnode)

# Convert absolute node numbers to internal-node indices (1…Nnode)
to_idx   <- map[,1] - length(t_tt$tip.label)
from_idx <- map[,2] - length(t_iq$tip.label)

# Transfer bootstrap values where available
valid <- which(from_idx >= 1 & from_idx <= length(iq_bs) & to_idx >= 1 & to_idx <= length(tt_bs))
tt_bs[to_idx[valid]] <- iq_bs[from_idx[valid]]

# Attach to timetree as node labels
t_tt$node.label <- as.character(round(tt_bs))

cat("Mapped bootstrap to timetree internal nodes:",
    sum(!is.na(tt_bs)), "of", length(tt_bs), "\n")

# Save a new timetree WITH bootstrap labels
getwd()
file.exists("treetime_M")
dir.create("treetime_M", showWarnings = FALSE)  # ensure folder exists

out_with_bs <- file.path(getwd(), "treetime_M", "timetree_with_bootstrap.nexus")
write.nexus(t_tt, file = out_with_bs)

cat("Wrote:", out_with_bs, "\n")


## Read input data
tree     <- read.nexus("treetime_M/timetree_with_bootstrap.nexus")
metadata <- read_csv("../metadata.filtered.csv", show_col_types = FALSE)

## Harmonize metadata and parse date
metadata <- metadata %>%
  mutate(
    Segment = toupper(Segment),
    date = ifelse(is.na(date), paste0(Year, "-06-15"), as.character(date)),
    date = as.Date(date)
  )

## Ensure a join key that matches tree tip labels
tip_labels <- tree$tip.label
if (!"name" %in% names(metadata)) {
  if ("taxa" %in% names(metadata)) {
    metadata <- metadata %>% mutate(name = taxa)
  } else if ("label" %in% names(metadata)) {
    metadata <- metadata %>% mutate(name = label)
  } else {
    metadata <- metadata %>% mutate(name = paste0(accession, "|", date))
  }
}

## Standardize Host labels and keep relevant columns
metadata <- metadata %>%
  mutate(
    Host = case_when(
      Host %in% c("Human", "Homo sapiens") ~ "Homo sapiens",
      Host == "Livestock" ~ "Livestock",
      Host == "Wildlife"  ~ "Wildlife",
      Host == "Vector"    ~ "Vector",
      TRUE ~ Host
    )
  ) %>%
  select(taxa, accession, date, Year, Country, Host, Segment, ProtMut, TotProtMut, Lineage)

## Filter to tips present in tree
metadata <- metadata %>% filter(taxa %in% tip_labels)

## Collapse countries into clear Regions only (no "Other")
metadata <- metadata %>%
  mutate(
    Region = case_when(
      Country %in% c("Kenya","Uganda","Tanzania","Somalia","Sudan","Rwanda","Burundi","Ethiopia") ~ "East Africa",
      Country %in% c("Mauritania","Senegal","Guinea","Burkina Faso","Niger")                      ~ "West Africa",
      Country %in% c("South Africa","Namibia","Zimbabwe","Madagascar","Angola")                   ~ "Southern Africa",
      Country %in% c("Central African Republic")                                                  ~ "Central Africa",
      Country %in% c("Egypt")                                                                     ~ "North Africa",
      TRUE ~ "Unknown"
    )
  )

## Determine most recent sampling date for calendar axis
mrsd <- tryCatch(max(metadata$date, na.rm = TRUE), error = function(e) NA)
if (is.na(mrsd)) mrsd <- as.Date("2020-06-15")

## Base time-scaled tree
p_base <- ggtree(tree, mrsd = mrsd) + theme_tree2()

## Host colors (distinct, circular points)
host_cols <- c(
  "Homo sapiens" = "#D55E00",  # orange
  "Livestock"    = "#0072B2",  # blue
  "Wildlife"     = "#009E73",  # green
  "Vector"       = "#CC79A7"   # magenta
)

# Read the NEW time-scaled tree that has bootstrap labels
tree_bs <- read.nexus(file.path(getwd(), "treetime_M", "timetree_with_bootstrap.nexus"))

# Recompute mrsd if needed
if (!exists("mrsd")) {
  mrsd <- tryCatch(max(metadata$date, na.rm = TRUE), error = function(e) NA)
  if (is.na(mrsd)) mrsd <- as.Date("2020-06-15")
}

# Now this line will work
p_host_base <- ggtree(tree_bs, mrsd = mrsd) + theme_tree2()

# Extract ggtree data to position support markers
node_df <- p_host_base$data
node_df$support <- suppressWarnings(as.numeric(node_df$label))

message("Internal nodes with numeric support in plotted tree: ",
        sum(!node_df$isTip & !is.na(node_df$support)))

# Host palette
host_cols <- c(
  "Homo sapiens" = "#D55E00",
  "Livestock"    = "#0072B2",
  "Wildlife"     = "#009E73",
  "Vector"       = "#CC79A7"
)

# Build plot with supports
p_host <- p_host_base %<+% metadata +
  geom_tree() +
  # Support dots (moderate)
  geom_point(
    data = subset(node_df, !isTip & !is.na(support) & support >= 70),
    aes(x = x, y = y), inherit.aes = FALSE,
    shape = 21, fill = "grey20", size = 1.6, stroke = 0
  ) +
  # Support numbers (strong)
  geom_text(
    data = subset(node_df, !isTip & !is.na(support) & support >= 90),
    aes(x = x + 0.0025, y = y, label = round(support)),
    inherit.aes = FALSE, size = 3.0, vjust = -0.15, color = "black"
  ) +
  geom_tippoint(aes(fill = Host), shape = 21, size = 5.6,
                color = "black", stroke = 0.35) +
  scale_fill_manual(values = host_cols, name = "Host") +
  scale_y_continuous(breaks = NULL, labels = NULL, expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(
    breaks = seq(1910, 2020, by = 10),
    labels = function(x) sprintf("%d", x),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Host-Associated Phylogenetic Structure of RVFV M Segment",
    x = "Year", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.9, "lines"),
    panel.grid = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    # >>> ADDED BELOW <<<
    axis.line.x  = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x  = element_text(size = 10, face = "bold")
  )


## Align metadata with tree order
tree_ordered_metadata <- metadata[match(tree$tip.label, metadata$taxa), ]

## Heatmap data: Region + Lineage
heatmap_data <- tree_ordered_metadata %>%
  mutate(y = 1:n()) %>%
  select(y, Region, Lineage)

## Region heatmap
## Region heatmap (custom, controlled colors)
p_region_hm <- ggplot(heatmap_data, aes(x = "Region", y = y, fill = Region)) +
  geom_tile() +
  scale_fill_manual(
    name = "Region",
    values = c(
      "East Africa"     = "#1B9E77",  # deep teal-green
      "West Africa"     = "#7570B3",  # muted violet
      "Southern Africa" = "#8C7853",  # dull brownish khaki
      "Central Africa"  = "#E6AB02",  # ochre yellow
      "North Africa"    = "#C2B280"   # desert beige / sand tone
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.8, "lines"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

## Lineage heatmap
p_lineage_hm <- ggplot(heatmap_data, aes(x = "Lineage", y = y, fill = Lineage)) +
  geom_tile() +
  scale_fill_manual(
    name = "Lineage",
    values = c(
      "A" = "#8B715A",  # warm taupe
      "B" = "#377EB8",  # blue
      "C" = "#4DAF4A",  # green
      "D" = "#984EA3",  # purple
      "E" = "#4169E1",  # royal blue
      "F" = "#FFFF33",  # yellow
      "G" = "#E97451",  # burnt sienna  
      "H" = "#0A4C6A", # dark cyan
      "I" = "#580F41",  # aubergine
      "K" = "#556B2F",  # deep olive green
      "L" = "#8DA0CB",   # lavender
      "N" = "#287C8E"   # teal blue
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_blank(),   # <--- added this line
    axis.ticks.y = element_blank()    # <--- added this line
  )


## Combine plots (Tree + Region + Lineage)
combined_plot <- p_host + p_region_hm + p_lineage_hm +
  plot_layout(widths = c(0.7, 0.15, 0.15))

print(combined_plot)

## Save final outputs
ggsave("M_combined_time_tree_regions.pdf",
       plot = combined_plot, width = 16, height = 10,
       units = "in", limitsize = FALSE)

ggsave("M_combined_time_tree_regions.png",
       plot = combined_plot, width = 16, height = 10,
       units = "in", dpi = 600, bg = "white")

## (Optional) Region × Lineage frequency summary
region_lineage_summary <- metadata %>%
  count(Region, Lineage) %>%
  arrange(Region, desc(n))
write_csv(region_lineage_summary, "Region_Lineage_summary.csv")
