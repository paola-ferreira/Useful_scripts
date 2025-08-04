# This script reads SRA run accession numbers from a text file ("runs_groups.txt"),
# fetches the read counts from NCBI for each run,
# sums the reads for groups of runs,
# and summarizes the associated BioProject accessions.

from Bio import Entrez
import re

Entrez.email = "your_email_here@alumni.usp.br"

def fetch_run_info(run_id):
    time.sleep(0.4)  # delay to avoid hitting NCBI too fast
    try:
        handle = Entrez.efetch(db="sra", id=run_id, rettype="xml", retmode="text")
        data = handle.read().decode("utf-8")  # decode bytes to string here
        handle.close()
    except Exception as e:
        raise RuntimeError(f"Entrez fetch failed for {run_id}: {e}")

    total_spots_match = re.search(r'total_spots="(\d+)"', data)
    total_spots = int(total_spots_match.group(1)) if total_spots_match else 0

    bioproject_match = re.search(r'namespace="BioProject">(\w+)', data)
    bioproject = bioproject_match.group(1) if bioproject_match else "Unknown"

    return total_spots, bioproject

def process_runs_group(run_ids):
    total_reads = 0
    bioprojects = set()

    for run_id in run_ids:
        if not run_id.startswith(('ERR', 'SRR', 'DRR')):
            print(f"[!] Skipping invalid run ID: {run_id}")
            continue
        try:
            reads, bioproject = fetch_run_info(run_id)
            print(f"  {run_id}: reads={reads}, BioProject={bioproject}")
            total_reads += reads
            bioprojects.add(bioproject)
        except Exception as e:
            print(f"[!] Error fetching {run_id}: {e}")

    return total_reads, bioprojects

def process_groups_from_file(input_file, output_file="groups_summary.txt"):
    with open(input_file) as f, open(output_file, "w") as out_f:
        for line_num, line in enumerate(f, 1):
            run_ids = line.strip().replace(',', ' ').split()
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

    print(f"All groups processed. Summary saved to {output_file}")

if __name__ == "__main__":
    process_groups_from_file("runs_groups.txt")

