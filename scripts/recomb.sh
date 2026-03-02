## Recombination detection was performed using the datamonkey platform: 
https://www.datamonkey.org/ 

#  Codon alignment for each of the RVFV genome fragments (S, M, L)

#!/usr/bin/env bash
set -euo pipefail

# Inputs
S_IN="S.fasta"
M_IN="M.fasta"
L_IN="L.fasta"

# Outputs
S_OUT="S.mafft.aln.fasta"
M_OUT="M.mafft.aln.fasta"
L_OUT="L.mafft.aln.fasta"

# MAFFT alignments
mafft --auto --reorder "${S_IN}" > "${S_OUT}"
mafft --auto --reorder "${M_IN}" > "${M_OUT}"
mafft --auto --reorder "${L_IN}" > "${L_OUT}"

echo "Done:"
echo "  ${S_OUT}"
echo "  ${M_OUT}"
echo "  ${L_OUT}"

chmod +x recomb.sh
./recomb.sh

# Each codon alignment was uploaded to the Datamonkey web server
  and analyzed using the Genetic Algorithm for Recombination Detection (GARD). 
  
# Breakpoint model comparison: Evaluated using ΔAIC-c
Model comparison based on corrected Akaike Information Criterion (ΔAIC-c) 
