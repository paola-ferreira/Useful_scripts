## This script roots all phylogenetic trees in a given directory using a specified outgroup taxon. 
## Input: unrooted Newick tree files (.treefile). 
## Output: rooted Newick tree files saved in the chosen output directory, preserving the original filenames
## Script written by Paola de L. Ferreira on Feb 14th, 2024

from Bio import Phylo
import os

# Specify the directory containing the input tree files
input_directory = "/home/paola/faststorage/17.Final_organization/6.Merging_subfamilies/1.Constrained/12.Rooting_genes_only_with_Dasypogon_included/"

# Specify the outgroup taxon
outgroup_taxon = "Dasypogon01"

# Specify the directory to save the rooted trees
output_directory = "/home/paola/faststorage/17.Final_organization/6.Merging_subfamilies/1.Constrained/13.Genes_trees_rooted_in_Dasypogon/"

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Loop through each file in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith(".treefile"):
        # Load the phylogenetic tree
        tree = Phylo.read(os.path.join(input_directory, filename), "newick")

        # Find the outgroup clade
        outgroup_clade = None
        for clade in tree.find_clades():
            if clade.name is not None and outgroup_taxon in clade.name:
                outgroup_clade = clade
                break

        # Root the tree using the outgroup
        if outgroup_clade:
            tree.root_with_outgroup(outgroup_clade)

        # Write the rooted tree to a new file in the output directory
        output_filename = os.path.join(output_directory, filename)
        Phylo.write(tree, output_filename, "newick")

        print(f"Rooted tree saved to: {output_filename}")
