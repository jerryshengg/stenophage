#!/usr/bin/env python3

"""
rscu.py

Uses FASTA files from:
~/Stenophage/Phages/fasta

Also includes:
~/Stenophage/RSCU/K279a.fasta
~/Stenophage/RSCU/D457.fasta

Uses GenBank files from:
~/Stenophage/Phages

Outputs to ~/Stenophage/RSCU:
codon_usage_counts.csv
codon_usage_RSCU.csv
trna_long_table.csv
trna_totals_by_phage.csv
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Data import CodonTable
import csv
import re
import math
import collections

fasta_folder = Path("~/Stenophage/Phages/fasta").expanduser()
gb_folder = Path("~/Stenophage/Phages").expanduser()
output_folder = Path("~/Stenophage/RSCU").expanduser()
output_folder.mkdir(parents=True, exist_ok=True)

CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
STD_TABLE = CodonTable.unambiguous_dna_by_id[11]

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V","SEC":"U","PYL":"O"
}

RE_TRNA = re.compile(r"tRNA[-\s]?([A-Za-z]{3}|[A-Z])(?:.*?\banticodon\b[:\s]*([ACGTUNWSMKRYBDHV]{3}))?", re.IGNORECASE)

def name_from_file(path):
    return path.stem.replace(" ", "_")

def trim_to_frame(seq):
    r = len(seq) % 3
    return seq if r == 0 else seq[:-r]

def count_codons(dna):
    dna = dna.upper().replace("U", "T")
    cnt = collections.Counter()
    for i in range(0, len(dna), 3):
        cod = dna[i:i + 3]
        if len(cod) == 3 and all(c in "ACGT" for c in cod):
            cnt[cod] += 1
    return cnt

def build_tables(table_id=11):
    tbl = CodonTable.unambiguous_dna_by_id.get(table_id, STD_TABLE)
    codon2aa = {}
    for cod in CODONS:
        if cod in tbl.stop_codons:
            codon2aa[cod] = "*"
        elif cod in tbl.forward_table:
            codon2aa[cod] = tbl.forward_table[cod]
        else:
            codon2aa[cod] = "*"

    aa2codons = collections.defaultdict(list)
    for cod, aa in codon2aa.items():
        if aa != "*":
            aa2codons[aa].append(cod)
    return codon2aa, aa2codons

CODON2AA, AA2CODONS = build_tables(11)

def gc3_from_counts(counts):
    total = sum(counts.values())
    if total == 0:
        return math.nan
    return sum(n for cod, n in counts.items() if cod[2] in "GC") / total

def rscu_from_counts(counts):
    rscu = {}
    for aa, codon_list in AA2CODONS.items():
        k = len(codon_list)
        total_for_aa = sum(counts.get(c, 0) for c in codon_list)

        if total_for_aa == 0:
            for c in codon_list:
                rscu[c] = math.nan
        else:
            expected = total_for_aa / k
            for c in codon_list:
                rscu[c] = counts.get(c, 0) / expected

    for c in CODONS:
        if c not in rscu:
            rscu[c] = math.nan

    return rscu

def parse_trna_feature(feat):
    aa_one = None
    anti = None
    prod = None

    if "product" in feat.qualifiers:
        prod = " ".join(feat.qualifiers["product"])
    elif "gene" in feat.qualifiers:
        prod = " ".join(feat.qualifiers["gene"])

    if prod:
        m = RE_TRNA.search(prod)
        if m:
            aa = m.group(1)
            if aa:
                aaU = aa.upper()
                if len(aaU) == 1:
                    aa_one = aaU
                elif len(aaU) == 3 and aaU in AA3_TO_1:
                    aa_one = AA3_TO_1[aaU]
            if m.group(2):
                anti = m.group(2).upper().replace("U", "T")

    if "anticodon" in feat.qualifiers and not anti:
        txt = " ".join(feat.qualifiers["anticodon"])
        m2 = re.search(r"seq:([ACGTUNWSMKRYBDHV]{3})", txt)
        if m2:
            anti = m2.group(1).upper().replace("U", "T")

        if not aa_one:
            m3 = re.search(r"aa:([A-Za-z]{3}|[A-Z])", txt)
            if m3:
                aaU = m3.group(1).upper()
                aa_one = aaU if len(aaU) == 1 else AA3_TO_1.get(aaU)

    return aa_one or "?", anti or "?"

def write_csv(path, fieldnames, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

# RSCU / codon usage from FASTA files

fasta_files = sorted(
    list(fasta_folder.glob("*.fasta")) +
    list(fasta_folder.glob("*.fa")) +
    list(fasta_folder.glob("*.fna"))
)

fasta_files = [f for f in fasta_files if f.name != "cog_phages.fasta"]

for extra in [output_folder / "K279a.fasta", output_folder / "D457.fasta"]:
    if extra.exists():
        fasta_files.append(extra)
        print("Added host genome:", extra.name)
    else:
        print("Missing host genome:", extra.name)

codon_counts_rows = []
rscu_rows = []

print(f"Found {len(fasta_files)} FASTA files")

for fasta in fasta_files:
    name = name_from_file(fasta)
    codon_counts = collections.Counter()

    print("Processing FASTA:", fasta.name)

    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq).upper().replace("U", "T")
        seq = trim_to_frame(seq)
        codon_counts += count_codons(seq)

    total_codons = sum(codon_counts.values())
    gc3 = gc3_from_counts(codon_counts)
    rscu = rscu_from_counts(codon_counts)

    counts_row = {"phage": name, "total_codons": total_codons, "GC3": gc3}
    for cod in CODONS:
        counts_row[cod] = codon_counts.get(cod, 0)
    codon_counts_rows.append(counts_row)

    rscu_row = {"phage": name}
    for cod in CODONS:
        rscu_row[cod] = rscu[cod]
    rscu_rows.append(rscu_row)

write_csv(output_folder / "codon_usage_counts.csv", ["phage", "total_codons", "GC3"] + CODONS, codon_counts_rows)
write_csv(output_folder / "codon_usage_RSCU.csv", ["phage"] + CODONS, rscu_rows)

# tRNA tables from GenBank files

gb_files = sorted(
    list(gb_folder.glob("*.gb")) +
    list(gb_folder.glob("*.gbk")) +
    list(gb_folder.glob("*.gbff"))
)

trna_long_rows = []
trna_totals_rows = []

print(f"Found {len(gb_files)} GenBank files")

for gb in gb_files:
    name = name_from_file(gb)
    trna_by_key = collections.Counter()
    trna_by_aa = collections.Counter()

    print("Processing GenBank:", gb.name)

    for rec in SeqIO.parse(gb, "genbank"):
        for feat in rec.features:
            if feat.type.lower() == "trna":
                aa, anti = parse_trna_feature(feat)
                trna_by_key[(aa, anti)] += 1
                trna_by_aa[aa] += 1

    for (aa, anti), n in sorted(trna_by_key.items()):
        trna_long_rows.append({
            "phage": name,
            "aa": aa,
            "anticodon": anti,
            "count": n
        })

    row = {
        "phage": name,
        "tRNA_total": sum(trna_by_key.values())
    }

    for aa in sorted(trna_by_aa):
        row[f"tRNA_{aa}"] = trna_by_aa[aa]

    trna_totals_rows.append(row)

aa_cols = sorted({k for r in trna_totals_rows for k in r if k.startswith("tRNA_") and k != "tRNA_total"})

write_csv(output_folder / "trna_long_table.csv", ["phage", "aa", "anticodon", "count"], trna_long_rows)
write_csv(output_folder / "trna_totals_by_phage.csv", ["phage", "tRNA_total"] + aa_cols, trna_totals_rows)

print("\nSaved:")
print(output_folder / "codon_usage_counts.csv")
print(output_folder / "codon_usage_RSCU.csv")
print(output_folder / "trna_long_table.csv")
print(output_folder / "trna_totals_by_phage.csv")
