#!/bin/bash
## This script checks FASTA alignment files to ensure they contain 
## at least 4 species (required for quartet-based analyses in ASTRAL).  
## For each file, it counts the number of unique species (based on FASTA headers) 
## and reports any alignments with fewer than 4 species.

# Loop through all .fasta files in the current directory
for fasta_file in *.fasta; do
    # Count the number of species in the alignment
    species_count=$(grep -c ">" "$fasta_file")

    # Check if the species count is less than 4 and print filename
    if [ "$species_count" -lt 4 ]; then
        echo "Alignment in $fasta_file has less than 4 species"
    fi
done
