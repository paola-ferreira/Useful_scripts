# This script reads SRA run accession numbers from a text file ("runs_groups.txt"),
# fetches the read counts from NCBI for each run,
# sums the reads for groups of runs,
# and summarizes the associated BioProject accessions.

from Bio import Entrez
import re

Entrez.email = "your_email_here@alumni.usp.br"

def fetch_run_info(run_id):
    handle = Entrez.efetch(db="sra", id=run_id, rettype="gb", retmode="text")
    data = handle.read()
    handle.close()
    if isinstance(data, bytes):
        data = data.decode('utf-8')

    total_spots_match = re.search(r'total_spots="(\d+)"', data)
    total_spots = int(total_spots_match.group(1)) if total_spots_match else 0

    bioproject_match = re.search(r'namespace="BioProject">(\w+)', data)
    bioproject = bioproject_match.group(1) if bioproject_match else "Unknown"

    return total_spots, bioproject

def process_runs_group(run_ids):
    total_reads = 0
    bioprojects = set()

    for run_id in run_ids:
        reads, bioproject = fetch_run_info(run_id)
        total_reads += reads
        bioprojects.add(bioproject)

    return total_reads, bioprojects

def main(input_file):
    with open(input_file) as f, open("groups_summary.txt", "w") as out_f:
        for line_num, line in enumerate(f, 1):
            run_ids = line.strip().split()
            if not run_ids:
                continue

            print(f"Processing group {line_num}: {run_ids}")
            total_reads, bioprojects = process_runs_group(run_ids)

            bioprojects_str = ", ".join(sorted(bioprojects))
            summary = (f"Group {line_num}: Runs: {', '.join(run_ids)}\n"
                       f"  Total reads: {total_reads}\n"
                       f"  BioProjects: {bioprojects_str}\n\n")

            print(summary)
            out_f.write(summary)

    print("All groups processed. Summary saved to groups_summary.txt")

if __name__ == "__main__":
    main("runs_groups.txt")

