# Script to remove alignments with fewer than 04 species

import os
from pathlib import Path

base_dir = Path(".")
low_species_file = "low_species.txt"
deleted_files = []

# Read all alignments to delete
with open(low_species_file) as f:
    alignments = [line.strip() for line in f]

for fasta_name in alignments:
    # Extract gene prefix up to _consensus
    prefix = fasta_name.split("_consensus")[0] + "_consensus"

    # Search directories: current + 0.* folders
    search_dirs = [base_dir] + [Path(d) for d in os.listdir(base_dir) if d.startswith("0.") and Path(d).is_dir()]

    for directory in search_dirs:
        for file_path in directory.glob(f"{prefix}*"):
            if file_path.is_file():
                try:
                    os.remove(file_path)
                    deleted_files.append(str(file_path.resolve()))
                except Exception as e:
                    print(f" Could not remove {file_path}: {e}")

print(" Cleanup complete.")
print(f"Total deleted files: {len(deleted_files)}")

# Optional: save log
if deleted_files:
    with open("deleted_files.log", "w") as log:
        log.write("\n".join(deleted_files))
