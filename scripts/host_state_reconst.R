
- Host-state reconstruction & spillover on RVFV M segment   
- Time-scaled tree (M segment)                            
- Ancestral host reconstruction (ER model)                
- Host→host transitions + network plot                    
###############################################################

# Load packages ---
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ape)
  library(phytools)
  library(igraph)
  library(scales)    # for rescale()
})

# User-editable paths (relative to repo root) ---
TREE_FILE <- "data/timetree_with_bootstrap.nexus"
META_FILE <- "data/metadata.csv"
OUTDIR    <- "outputs/host_state_M"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Input checks ---
stopifnot(file.exists(TREE_FILE))
stopifnot(file.exists(META_FILE))

# Load time-scaled tree & metadata ---
tree_bs  <- read.nexus(TREE_FILE)
metadata <- read_csv(META_FILE, show_col_types = FALSE)

required_cols <- c("taxa", "Host", "Segment")
missing_cols  <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("metadata.csv is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Clean metadata: standardise Host, keep only M segment ---
metadata_M <- metadata %>%
  mutate(
    Segment = toupper(Segment),
    Host = case_when(
      Host %in% c("Human", "Homo sapiens", "human") ~ "human",
      Host %in% c("Livestock", "livestock")         ~ "livestock",
      Host %in% c("Wildlife", "wildlife")           ~ "wildlife",
      Host %in% c("Vector", "vector")               ~ "vector",
      TRUE                                          ~ NA_character_
    )
  ) %>%
  filter(
    Segment == "M",
    Host %in% c("human", "livestock", "wildlife", "vector"),
    !is.na(taxa)
  )

if (nrow(metadata_M) == 0) {
  stop("No rows remain after filtering for Segment == 'M' and valid Host categories.")
}

# Prepare tree for ACE: root & resolve polytomies ---
tree_fix <- tree_bs %>%
  midpoint.root() %>%
  multi2di()

# Prune tree to tips that have metadata ---
tips_to_keep <- intersect(tree_fix$tip.label, metadata_M$taxa)
if (length(tips_to_keep) < 5) {
  stop("Too few matched tips between tree and metadata (n = ", length(tips_to_keep), "). Check taxa names.")
}

tree_fix <- drop.tip(tree_fix, setdiff(tree_fix$tip.label, tips_to_keep))

# Reorder metadata to match tip order
metadata_ord <- metadata_M[match(tree_fix$tip.label, metadata_M$taxa), ]
if (any(is.na(metadata_ord$taxa))) {
  stop("Metadata reordering failed: some pruned tree tips have no matching metadata row.")
}

# Build tip state vector (host categories) ---
tip_states <- factor(metadata_ord$Host)
tip_states <- droplevels(tip_states)
names(tip_states) <- tree_fix$tip.label

# Ancestral host reconstruction (equal-rates ER model) ---
host_ace <- ace(tip_states, tree_fix, type = "discrete", model = "ER")

# Get most likely host for each internal node ---
internal_states <- apply(host_ace$lik.anc, 1, function(x) names(which.max(x)))
ntips <- length(tree_fix$tip.label)

get_host_state <- function(node) {
  if (node <= ntips) {
    as.character(tip_states[node])      # tip node
  } else {
    internal_states[node - ntips]       # internal node
  }
}

# Build host→host transition table ---
edge_df <- as.data.frame(tree_fix$edge)
colnames(edge_df) <- c("parent", "child")

edge_df <- edge_df %>%
  mutate(
    from_host = sapply(parent, get_host_state),
    to_host   = sapply(child,  get_host_state)
  )

host_transition_table <- edge_df %>%
  count(from_host, to_host, name = "n") %>%
  arrange(desc(n))

spillover_into_human <- host_transition_table %>%
  filter(to_host == "human", from_host != "human", n > 0) %>%
  arrange(desc(n))

# Save transition tables ---
write_csv(host_transition_table, file.path(OUTDIR, "Host_Transitions_M_segment.csv"))
write_csv(spillover_into_human, file.path(OUTDIR, "Spillover_to_Humans_M_segment.csv"))

# Plot host-to-host transition network ---
# Ensure fixed vertex ordering (even if some states absent)
vertices <- data.frame(
  name  = c("human", "livestock", "wildlife", "vector"),
  color = c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
)

g <- graph_from_data_frame(
  host_transition_table %>% rename(from = from_host, to = to_host),
  directed = TRUE,
  vertices = vertices
)

layout_circle <- layout_in_circle(g) * 0.8

# Edge aesthetics
E(g)$width      <- rescale(E(g)$n, to = c(1, 4))
E(g)$arrow.size <- 0.5
E(g)$color      <- "grey30"

# Node aesthetics
V(g)$size        <- 40
V(g)$label       <- V(g)$name
V(g)$label.cex   <- 1
V(g)$label.color <- "black"
V(g)$color       <- V(g)$color

png(file.path(OUTDIR, "Host_Transition_Network.png"), width = 2000, height = 2000, res = 300)
plot(
  g,
  layout      = layout_circle,
  edge.curved = 0.2,
  margin      = 0.3,
  main        = "Host-to-Host Transition Network (RVFV M Segment)"
)
dev.off()

# Save session info (reproducibility) ---
sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. Outputs written to: ", OUTDIR)
                         

