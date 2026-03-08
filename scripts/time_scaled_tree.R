### Phylogentic Inference:

RVFV M segment 
# Tools needed:
- MAFFT    # Multiple Sequence Alignment
- trimAl   # trims poorly aligned regions
- IQ-TREE  # Maximum likelihood phylogenetic tree inference
- TreeTime # time scaled phylogenetic reconstruction

## R part
# Set the WD
setwd("~/IsaacOmara/Desktop/Manuscripts/M") #Update accordingly

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


#### RVFV L Segment

# Get WD
getwd()

# Set the WD
setwd("~/IsaacOmara/Desktop/Manuscripts/L")

# Load the libraries
library(readr)
library(dplyr)
library(tibble)

# Read inputs 
human  <- read_csv("sequences.csv",  show_col_types = FALSE)
animal <- read_csv("sequences1.csv", show_col_types = FALSE)

# Merge, filter L, standardize Host, build date & taxa --------
meta_L <- bind_rows(human, animal) %>%
  filter(toupper(Segment) == "L") %>%
  mutate(
    Host = case_when(
      Host %in% c("Human","Homo sapiens") ~ "Homo sapiens",
      Host == "Livestock" ~ "Livestock",
      Host == "Wildlife"  ~ "Wildlife",
      Host == "Vector"    ~ "Vector",
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

# Add a single L reference 
ref_acc_L <- "NC_014397"
ref_date  <- as.Date("1977-06-15")

ref_row <- tibble(
  taxa       = paste0(ref_acc_L, "|", ref_date),
  accession  = ref_acc_L,
  date       = ref_date,
  Year       = as.integer(format(ref_date, "%Y")),
  Country    = "Reference",
  Host       = "Reference",
  Segment    = "L",
  ProtMut    = NA_character_,
  TotProtMut = NA_real_,
  Lineage    = "A"
)

meta_L <- bind_rows(meta_L, ref_row) %>% distinct(taxa, .keep_all = TRUE)

# Write metadata + TreeTime helpers --------
write_csv(meta_L, "metadata.csv")

write.table(meta_L %>% select(name = taxa, date),
            file = "dates_L.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

writeLines(meta_L$taxa, "ids_L.txt")

cat("L metadata rows:", nrow(meta_L), "\n")
print(head(meta_L, 3))

### BASH SCRIPTING: Relabel Fasta Headers with awk;
awk -F, 'NR==FNR{map[$2]=$1; next} /^>/{hdr=$0; sub(/^>/,"",hdr); split(hdr,a,/[\t |]/); acc=a[1];
  if(acc in map){print ">"map[acc]} else {print ">"hdr} next} {print}' \
metadata.csv sequences.fasta > human.relabel.fasta

awk -F, 'NR==FNR{map[$2]=$1; next} /^>/{hdr=$0; sub(/^>/,"",hdr); split(hdr,a,/[\t |]/); acc=a[1];
  if(acc in map){print ">"map[acc]} else {print ">"hdr} next} {print}' \
metadata.csv sequences1.fasta > nonhuman.relabel.fasta


## COMBINE, dedup, align & trim
seqkit grep -f ids_L.txt human.relabel.fasta nonhuman.relabel.fasta > riftL_all.fasta
seqkit rmdup -s riftL_all.fasta > riftL_all.nodup.fasta
grep '^>' riftL_all.nodup.fasta | sed 's/^>//' > kept_ids.txt

mafft --auto --reorder --adjustdirectionaccurately --thread -1 riftL_all.nodup.fasta > riftL_all_align.fasta
trimal -automated1 -in riftL_all_align.fasta -out riftL_all_align.trim.fasta

# Reference present?
grep '^>' riftL_all_align.fasta | grep '^>NC_014397|'
grep '^>' riftL_all_align.trim.fasta | grep '^>NC_014397|'

### IQ-Tree:
iqtree2 -s riftL_all_align.trim.fasta \
-m MFP -B 1000 -alrt 1000 -nt AUTO \
-o "NC_014397|1977-06-15" \
-pre L_all

## Tree-time/timetree.nexus
treetime --tree L_all.treefile \
--aln  riftL_all_align.trim.fasta \
--dates dates_L.tsv \
--clock-filter 3 \
--outdir treetime_L \
--gtr JC69

## UPDATE Metadata 

meta <- read_csv("metadata.csv", show_col_types = FALSE)
kept <- readr::read_lines("kept_ids.txt")
meta_filtered <- meta %>% filter(taxa %in% kept)
write_csv(meta_filtered, "metadata.filtered.csv")
cat("Final L metadata rows:", nrow(meta_filtered), "\n")


                                                            #### RVFV S segment 

setwd("~/IsaacOmara/Desktop/Manuscripts/S")

library(readr)
library(dplyr)
library(tibble)

# Read inputs 
human  <- read_csv("sequences.csv",  show_col_types = FALSE)
animal <- read_csv("sequences1.csv", show_col_types = FALSE)

# Merge, filter S, standardize Host, build date & taxa 
meta_S <- bind_rows(human, animal) %>%
  filter(toupper(Segment) == "S") %>%
  mutate(
    Host = case_when(
      Host %in% c("Human","Homo sapiens") ~ "Homo sapiens",
      Host == "Livestock" ~ "Livestock",
      Host == "Wildlife"  ~ "Wildlife",
      Host == "Vector"    ~ "Vector",
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

# Add a single S reference (edit date if you know the true reference date) --------
ref_row <- tibble(
  taxa       = "NC_014395|1977-06-15",
  accession  = "NC_014395",
  date       = as.Date("1977-06-15"),
  Year       = 1977L,
  Country    = "Reference",
  Host       = "Reference",   # set to NA if you don't want it in legends
  Segment    = "S",
  ProtMut    = NA_character_,
  TotProtMut = NA_real_,
  Lineage    = "A"
)

meta_S <- bind_rows(meta_S, ref_row) %>% distinct(taxa, .keep_all = TRUE)

# Write metadata + TreeTime helpers --------
write_csv(meta_S, "metadata.csv")

write.table(meta_S %>% select(name = taxa, date),
            file = "dates_S.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

writeLines(meta_S$taxa, "ids_S.txt")

cat("S metadata rows:", nrow(meta_S), "\n")
print(head(meta_S, 3))

## RENAMING HEADERS TO MATCH Taxa: Bash Scripting
# Build mapping dictionary in awk and substitute headers
awk -F, 'NR==FNR{map[$2]=$1; next} /^>/{acc=$0; sub(/^>/,"",acc); if(acc in map){print ">"map[acc]} else {print ">"acc}} !/^>/' \
metadata.csv sequences.fasta > human.relabel.fasta

awk -F, 'NR==FNR{map[$2]=$1; next} /^>/{acc=$0; sub(/^>/,"",acc); if(acc in map){print ">"map[acc]} else {print ">"acc}} !/^>/' \
metadata.csv sequences1.fasta > nonhuman.relabel.fasta

# Cross check
grep "^>" human.relabel.fasta | head
grep "^>" nonhuman.relabel.fasta | head

# And should look like
>MN123456|2016-06-15
>AB987654|2007-06-15

# Build the Combined Fasta for both hosts using ids_L.txt
seqkit grep -f ids_S.txt human.relabel.fasta nonhuman.relabel.fasta > riftS_all.fasta

# QC: remove exact duplicates, filter extreme lengths if needed
seqkit rmdup -s riftS_all.fasta > riftS_all.nodup.fasta # [83 duplicates removed =146]

# capture the IDs that survived
grep '^>' riftS_all.nodup.fasta | sed 's/^>//' > kept_ids.txt

# tweak reference length / tolerance to your M segment if you want strict filtering
mafft --auto --reorder --adjustdirectionaccurately --thread -1 riftS_all.nodup.fasta > riftS_all_align.fasta
trimal -automated1 -in riftS_all_align.fasta -out riftS_all_align.trim.fasta

# Check whether the reference is still there; 
# After MAFFT
grep '^>' riftS_all_align.fasta | grep '^>NC_014395|'

# After trimAl
grep '^>' riftS_all_align.trim.fasta | grep '^>NC_014395|'

## ML Tree -IQ-Tree
iqtree2 -s riftS_all_align.trim.fasta \
-m MFP -B 1000 -alrt 1000 -nt AUTO \
-o "NC_014395|1977-06-15" \
-pre S_all
# S_all.treefile will be written rooted on that taxon

#### TIME-SCALE TREE = treetime.nexus
treetime --tree S_all.treefile \
--aln  riftS_all_align.trim.fasta \
--dates dates_S.tsv \
--clock-filter 3 \
--outdir treetime_S \
--gtr JC69
# treetime_S/timetree.nexus


#### UPDATING Metadata to match fasta sequence remaining after removing duplicates
meta <- read_csv("metadata.csv") # 228
kept <- read_lines("kept_ids.txt")

meta_filtered <- meta %>% filter(taxa %in% kept) # 146
write_csv(meta_filtered, "metadata.filtered.csv") 
