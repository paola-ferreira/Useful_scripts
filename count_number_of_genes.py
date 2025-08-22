#!/bin/bash
# This script iterates through all FASTA files in a specified directory, 
# counts the number of gene entries (identified by header lines starting with '>') 
# in each file, and outputs the results as a summary table (species + gene count) 
# into a text file.

# Define the directory containing your FASTA files
fasta_directory="/home/paola/faststorage/17.Final_organization/3.Arecoids/1.Constrained/5.Coverage"

# Define the output file
output_file="gene_counts.txt"

# Remove existing output file if it exists
rm -f "$output_file"

# Loop through each FASTA file in the directory
for fasta_file in "$fasta_directory"/*.fasta; do
    if [ -f "$fasta_file" ]; then
        # Get the file name without extension
        file_name=$(basename "$fasta_file" .fasta)

        # Count the number of genes and append the result to the output file
        gene_count=$(sed -n '/^>/p' "$fasta_file" | wc -l)
        echo "$file_name $gene_count" >> "$output_file"
    fi
done

echo "Gene counts saved to $output_file"
