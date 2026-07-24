#!/usr/bin/env python3

"""
compile_phages.py

Combines all FASTA files from:
~/Stenophage/Phages/fasta/

Output:
~/Stenophage/Phages/fasta/cog_phages.fasta
"""

from pathlib import Path
from Bio import SeqIO

input_folder = Path("~/Stenophage/Phages/fasta").expanduser()
output_fasta = input_folder / "cog_phages.fasta"

fasta_files = sorted(
    list(input_folder.glob("*.fasta")) +
    list(input_folder.glob("*.fa")) +
    list(input_folder.glob("*.fna"))
)

records = []

print(f"Looking in: {input_folder}")
print(f"Found {len(fasta_files)} FASTA files")

for fasta in fasta_files:

    if fasta.name == output_fasta.name:
        continue

    print("Adding:", fasta.name)

    for rec in SeqIO.parse(fasta, "fasta"):
        records.append(rec)

SeqIO.write(records, output_fasta, "fasta")

print("\nSaved:", output_fasta)
print("Sequences written:", len(records))
