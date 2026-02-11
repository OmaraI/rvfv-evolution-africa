### Phylogentic Inference/Reconstruction:RVFV M segment 
# Tools needed:
- MAFFT  
- trimAl
- IQ-TREE
- TreeTime

## R part
# Set the Working directory
setwd("/Users/IsaacOmara/Desktop/Manuscripts/Data/R-Work/Monika/Nexus-Files/M")

# Load libraries
library(readr)
library(dplyr)
library(tibble)

# Read inputs 
human  <- read_csv("sequences.csv",  show_col_types = FALSE)
animal <- read_csv("sequences1.csv", show_col_types = FALSE)

# Merge, filter M, standardize Host build date & taxa 
meta_M <- bind_rows(human, animal) %>%
  filter(toupper(Segment) == "M") %>%
  mutate(
    Host = case_when(
      Host %in% c("Human","Homo sapiens") ~ "Homo sapiens",
      Host == "Livestock" ~ "Livestock",
      Host == "Wildlife"  ~ "Wildlife",
      Host == "Vector"    ~ "Vector"
    ),
    Year = as.integer(`Year of Sample Collection`),
    date = as.Date(paste0(Year, "-06-15")),
    taxa = paste0(accession, "|", date)
  ) %>%
  transmute(
    taxa, accession, date, Year, Country, Host, Segment,
    ProtMut    = `Protein Mutations`,
    TotProtMut = `Total Protein Mutations`,
    Lineage
  ) %>%
  filter(!is.na(accession), !is.na(date)) %>%
  distinct(taxa, .keep_all = TRUE)

# Add a single M reference 
ref_acc_M <- "NC_014396"
ref_date  <- as.Date("1977-06-15")

ref_row <- tibble(
  taxa       = paste0(ref_acc_M, "|", ref_date),
  accession  = ref_acc_M,
  date       = ref_date,
  Year       = as.integer(format(ref_date, "%Y")),
  Country    = "Reference",
  Host       = "Reference",   # set to NA if you don't want it in legends
  Segment    = "M",
  ProtMut    = NA_character_,
  TotProtMut = NA_real_,
  Lineage    = "A"
)

meta_M <- bind_rows(meta_M, ref_row) %>% distinct(taxa, .keep_all = TRUE)

# Write metadata + TreeTime helpers --------
write_csv(meta_M, "metadata.csv")

write.table(meta_M %>% select(name = taxa, date),
            file = "dates_M.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

writeLines(meta_M$taxa, "ids_M.txt")

cat("M metadata rows:", nrow(meta_M), "\n")
print(head(meta_M, 3))

## BASH SCRIPTING 
# Bash Scripting: Relabel FASTA Headers using awk; 
awk -F, 'NR==FNR{map[$2]=$1; next} /^>/{hdr=$0; sub(/^>/,"",hdr); split(hdr,a,/[\t |]/); acc=a[1];
  if(acc in map){print ">"map[acc]} else {print ">"hdr} next} {print}' \
metadata.csv sequences.fasta > human.relabel.fasta

awk -F, 'NR==FNR{map[$2]=$1; next} /^>/{hdr=$0; sub(/^>/,"",hdr); split(hdr,a,/[\t |]/); acc=a[1];
  if(acc in map){print ">"map[acc]} else {print ">"hdr} next} {print}' \
metadata.csv sequences1.fasta > nonhuman.relabel.fasta

## Combine, dedup, align & trim
seqkit grep -f ids_M.txt human.relabel.fasta nonhuman.relabel.fasta > riftM_all.fasta
seqkit rmdup -s riftM_all.fasta > riftM_all.nodup.fasta
grep '^>' riftM_all.nodup.fasta | sed 's/^>//' > kept_ids.txt

## Multiple Sequence Alignment
mafft --auto --reorder --adjustdirectionaccurately --thread -1 riftM_all.nodup.fasta > riftM_all_align.fasta
trimal -automated1 -in riftM_all_align.fasta -out riftM_all_align.trim.fasta

# Reference present?
grep '^>' riftM_all_align.fasta | grep '^>NC_014396|'
grep '^>' riftM_all_align.trim.fasta | grep '^>NC_014396|'

### IQ-Tree + Tree time
iqtree2 -s riftM_all_align.trim.fasta \
-m MFP -B 1000 -bnni -alrt 1000 -nt AUTO \
-pre M_all

# Treetime/timetree.nexus
treetime --tree M_all.treefile \
--aln  riftM_all_align.trim.fasta \
--dates dates_M.tsv \
--clock-filter 3 \
--outdir treetime_M \
--gtr JC69

## Use R
### UPDATE Metadata to retain tips
meta <- read_csv("metadata.csv", show_col_types = FALSE)

kept <- readr::read_lines("kept_ids.txt")

meta_filtered <- meta %>% filter(taxa %in% kept)

write_csv(meta_filtered, "metadata.filtered.csv")

cat("Final M metadata rows:", nrow(meta_filtered), "\n")
