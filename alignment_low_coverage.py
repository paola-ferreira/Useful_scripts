## Script to identify alignments with fewer than 04 species

from pathlib import Path

base_dir = Path(".")
low_species_file = "low_species.txt"

with open(low_species_file, "w") as out_f:
    for fasta_path in base_dir.glob("*.fasta"):
        # Count species (lines starting with ">")
        species_count = sum(1 for line in open(fasta_path) if line.startswith(">"))

        if species_count < 4:
            print(f"{fasta_path.name} has {species_count} species")
            out_f.write(f"{fasta_path.name}\n")

print(f" Identification complete. See {low_species_file}")
