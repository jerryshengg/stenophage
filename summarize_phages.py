#!/usr/bin/env python3

"""
summarize_phages.py

Reads GenBank files from ~/Stenophage/Phages and extracts:
genome length, GC%, ORF/CDS counts, tRNA count, CDS density,
and basic CDS annotation classes.

Output:
~/Stenophage/Summary/phage_summary.csv
"""

from pathlib import Path
from Bio import SeqIO
import pandas as pd

input_folder = Path("~/Stenophage/Phages").expanduser()
output_csv = Path("~/Stenophage/Summary/phage_summary.csv").expanduser()

def gc_percent(seq):
    seq = str(seq).upper()
    if len(seq) == 0:
        return 0
    return ((seq.count("G") + seq.count("C")) / len(seq)) * 100

def is_pseudogene(feat):
    if "pseudo" in feat.qualifiers:
        return True
    if feat.qualifiers.get("pseudogene"):
        return True
    notes = feat.qualifiers.get("note", [])
    return any("pseudo" in n.lower() for n in notes)

def classify_cds(feat):
    text = " ".join(
        feat.qualifiers.get("product", []) +
        feat.qualifiers.get("function", []) +
        feat.qualifiers.get("note", [])
    ).lower()

    if "hypothetical" in text:
        return "hypothetical"
    if "putative" in text or "probable" in text or "predicted" in text:
        return "putative"
    return "functional"

rows = []

gb_files = sorted(input_folder.glob("*.gb*"))

print(f"Looking in: {input_folder}")
print(f"Found {len(gb_files)} files")

for file in gb_files:
    print("Processing:", file.name)

    for rec in SeqIO.parse(file, "genbank"):
        orf_count_total = 0
        cds_count = 0
        trna_count = 0
        cds_functional = 0
        cds_putative = 0
        cds_hypothetical = 0

        for feat in rec.features:
            if feat.type == "CDS":
                orf_count_total += 1

                if not is_pseudogene(feat):
                    cds_count += 1
                    cls = classify_cds(feat)

                    if cls == "functional":
                        cds_functional += 1
                    elif cls == "putative":
                        cds_putative += 1
                    elif cls == "hypothetical":
                        cds_hypothetical += 1

            elif feat.type.lower() == "trna":
                trna_count += 1

        rows.append({
            "phage": rec.id,
            "file": file.name,
            "genome_length_bp": len(rec.seq),
            "GC_percent": round(gc_percent(rec.seq), 2),
            "ORF_count_total": orf_count_total,
            "CDS_count": cds_count,
            "CDS_functional": cds_functional,
            "CDS_putative": cds_putative,
            "CDS_hypothetical": cds_hypothetical,
            "tRNA_count": trna_count,
            "CDS_per_10kb": round((cds_count / len(rec.seq)) * 10000, 2)
        })

df = pd.DataFrame(rows)

output_csv.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(output_csv, index=False)

print("\nSaved:", output_csv)
print("Rows:", len(df))

if len(df) > 0:
    print("Mean genome length:", round(df["genome_length_bp"].mean(), 1))
    print("Mean GC:", round(df["GC_percent"].mean(), 2))
    print("Mean total ORFs:", round(df["ORF_count_total"].mean(), 2))
    print("Mean CDS:", round(df["CDS_count"].mean(), 2))
    print("Mean tRNA:", round(df["tRNA_count"].mean(), 2))
else:
    print("No records parsed")