


## Host-state reconstruction & spillover analysis on RVFV M segment   
## - Time-scaled tree (M segment)                            
## - Ancestral host reconstruction (ER model)                
## - Host→host transitions + network plot                    


## 0. Load packages 
library(readr)
library(dplyr)
library(ape)
library(phytools)
library(igraph)
library(scales)

## 1. Set working directory 
setwd("~/Desktop/Manuscripts/Data/R-Work/Nexus-Files/M/Bootstrap/Host_State_R_SA/")

## 2. Load time-scaled tree & metadata 
tree_bs  <- read.nexus("timetree_with_bootstrap.nexus")
metadata <- read_csv("metadata.csv", show_col_types = FALSE)

## 3. Clean metadata: standardise Host, keep only M segment 
metadata_M <- metadata %>%
  mutate(
    Segment = toupper(Segment),
    Host = case_when(
      Host %in% c("Human", "Homo sapiens", "human") ~ "human",
      Host %in% c("Livestock", "livestock")         ~ "livestock",
      Host %in% c("Wildlife", "wildlife")           ~ "wildlife",
      Host %in% c("Vector", "vector")               ~ "vector",
      TRUE                                           ~ Host
    )
  ) %>%
  filter(
    Segment == "M",
    Host %in% c("human", "livestock", "wildlife", "vector"),
    !is.na(taxa)
  )

## 4. Prepare tree for ACE: root & resolve polytomies 
tree_fix <- tree_bs %>%
  midpoint.root() %>%
  multi2di()

## 5. Prune tree to tips that have metadata 
tips_to_keep <- intersect(tree_fix$tip.label, metadata_M$taxa)
tree_fix <- drop.tip(tree_fix, setdiff(tree_fix$tip.label, tips_to_keep))

metadata_ord <- metadata_M[match(tree_fix$tip.label, metadata_M$taxa), ]

## 6. Build tip state vector 
tip_states <- factor(metadata_ord$Host)
tip_states <- droplevels(tip_states)
names(tip_states) <- tree_fix$tip.label

## 7. Ancestral host reconstruction 
host_ace <- ace(tip_states, tree_fix, type = "discrete", model = "ER")

## 8. Get most likely host for each internal node 
internal_states <- apply(host_ace$lik.anc, 1, function(x) names(which.max(x)))
ntips <- length(tree_fix$tip.label)

get_host_state <- function(node) {
  if (node <= ntips) {
    as.character(tip_states[node])
  } else {
    internal_states[node - ntips]
  }
}

## 9. Build host→host transition table 
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

print(host_transition_table)

## 10. Extract non-human → human spillovers 
spillover_into_human <- host_transition_table %>%
  filter(to_host == "human", from_host != "human", n > 0) %>%
  arrange(desc(n))

print(spillover_into_human)

## 11. Save transition tables 
write_csv(host_transition_table, "Host_Transitions_M_segment.csv")
write_csv(spillover_into_human, "Spillover_to_Humans_M_segment.csv")


## 12. Plot host-to-host transition network

g <- graph_from_data_frame(
  host_transition_table,
  directed = TRUE,
  vertices = data.frame(
    name  = c("human", "livestock", "wildlife", "vector"),
    color = c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
  )
)

layout_circle <- layout_in_circle(g)
layout_circle <- layout_circle * 0.8

# General edge styling
E(g)$width       <- scales::rescale(as.numeric(E(g)$n), to = c(0.8, 2.5))
E(g)$arrow.size  <- 0.4
E(g)$color       <- "grey65"

# Fix self-loops (this is the key part)
loops <- which(is.loop(g))

E(g)$color[loops] <- "grey75"   # lighter than other edges
E(g)$width[loops] <- 0.8        # thinner loops
E(g)$label       <- as.character(E(g)$n)
E(g)$label.cex   <- 0.8
E(g)$label.color <- "black"

V(g)$size        <- 40
V(g)$label       <- V(g)$name
V(g)$label.cex   <- 1
V(g)$label.color <- "black"
V(g)$color       <- V(g)$color

## Plot on screen 
plot(
  g,
  layout             = layout_circle,
  edge.curved        = 0.25,
  edge.width         = E(g)$width,
  edge.arrow.size    = E(g)$arrow.size,
  edge.color         = E(g)$color,
  edge.label         = E(g)$label,
  edge.label.cex     = E(g)$label.cex,
  edge.label.color   = E(g)$label.color,
  edge.label.dist    = 0.6,
  edge.label.degree  = pi / 2,
  vertex.size        = V(g)$size,
  vertex.label       = V(g)$label,
  vertex.label.cex   = V(g)$label.cex,
  vertex.label.color = V(g)$label.color,
  vertex.color       = V(g)$color,
  vertex.frame.color = NA,   # emoves black borders
  margin             = 0.35
)

## Save high-resolution image 
png("Host_Transition_Network_with_n_labels.jpeg",
    width = 2000, height = 2000, res = 300)

plot(
  g,
  layout             = layout_circle,
  edge.curved        = 0.25,
  edge.width         = E(g)$width,
  edge.arrow.size    = E(g)$arrow.size,
  edge.color         = E(g)$color,
  edge.label         = E(g)$label,
  edge.label.cex     = E(g)$label.cex,
  edge.label.color   = E(g)$label.color,
  edge.label.dist    = 0.6,
  edge.label.degree  = pi / 2,
  vertex.size        = V(g)$size,
  vertex.label       = V(g)$label,
  vertex.label.cex   = V(g)$label.cex,
  vertex.label.color = V(g)$label.color,
  vertex.color       = V(g)$color,
  vertex.frame.color = NA,   
  margin             = 0.35
)

dev.off()

                         
