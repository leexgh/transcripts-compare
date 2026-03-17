#!/usr/bin/env python3
"""
Check RefSeq NM IDs from:
  - oncokb_isoform_versioned.tsv   (GRCh37 + GRCh38 columns)
  - isoform_overrides_at_mskcc_grch37.txt  (GRCh37 only)
  - Iv7_dmp_isoform_merged_overrides.txt   (GRCh37 only)

Against GFF files:
  - GCF_000001405.25_GRCh37.p13_genomic.gff  (GRCh37)
  - ref_GRCh37.p10_top_level.gff3            (GRCh37)
  - GCF_000001405.40_GRCh38.p14_genomic.gff  (GRCh38)

Logic:
  - "found" = exact versioned ID present  OR  base ID (NM_XXXXXX) present with any version
  - Reports IDs where NEITHER the exact version NOR any version of the base is found in any applicable GFF
"""
import re
import os
import sys
import gzip

ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(ROOT, "files")

def fp(*parts):
    return os.path.join(FILES, *parts)

NM_RE = re.compile(r"Name=(NM_\d+\.\d+)")

def _open(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")

# ---------------------------------------------------------------------------
# Build NM sets from a GFF file: {exact_id}, {base_id}
# ---------------------------------------------------------------------------
def build_nm_sets(gff_path: str, label: str):
    exact = set()
    bases = set()
    print(f"  Scanning {label}…", flush=True)
    with _open(gff_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            for m in NM_RE.finditer(line):
                nm = m.group(1)
                exact.add(nm)
                bases.add(nm.split(".")[0])
    print(f"    → {len(exact):,} versioned NMs, {len(bases):,} base IDs")
    return exact, bases

# ---------------------------------------------------------------------------
# Collect RefSeq IDs to check from source files
# ---------------------------------------------------------------------------
def collect_ids():
    """Returns list of (source_label, gene, nm_id, assembly)"""
    records = []

    # 1. oncokb_isoform_versioned.tsv — both assemblies
    tsv = fp("oncokb_isoform_versioned.tsv")
    with open(tsv, encoding="utf-8") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        def ci(name):
            for i, c in enumerate(cols):
                if c.strip() == name: return i
            return None
        ig = ci("Hugo Symbol"); i37r = ci("GRCh37 RefSeq"); i38r = ci("GRCh38 RefSeq")
        for line in fh:
            row = line.rstrip("\n").split("\t")
            gene = row[ig].strip() if ig is not None and ig < len(row) else ""
            if not gene: continue
            nm37 = row[i37r].strip() if i37r is not None and i37r < len(row) else ""
            nm38 = row[i38r].strip() if i38r is not None and i38r < len(row) else ""
            if nm37 and nm37.startswith("NM_"):
                records.append(("oncokb_isoform_versioned (GRCh37)", gene, nm37, "GRCh37"))
            if nm38 and nm38.startswith("NM_"):
                records.append(("oncokb_isoform_versioned (GRCh38)", gene, nm38, "GRCh38"))

    # 2. isoform_overrides_at_mskcc_grch37.txt — GRCh37
    # Format: gene_name  refseq_id  enst_id  note
    mskcc = fp("isoform_overrides_at_mskcc_grch37.txt")
    with open(mskcc, encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"): continue
            parts = line.split("\t")
            if len(parts) < 2: continue
            gene = parts[0].strip()
            nm   = parts[1].strip()
            if nm.startswith("NM_"):
                records.append(("mskcc_grch37", gene, nm, "GRCh37"))

    # 3. Iv7_dmp_isoform_merged_overrides.txt — GRCh37
    # Format: nm_id  gene  nm_base
    iv7 = fp("Iv7_dmp_isoform_merged_overrides.txt")
    with open(iv7, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"): continue
            parts = line.split("\t")
            if len(parts) < 1: continue
            nm = parts[0].strip()
            gene = parts[1].strip() if len(parts) > 1 else ""
            if nm.startswith("NM_"):
                records.append(("Iv7_overrides", gene, nm, "GRCh37"))

    return records

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=== Building NM sets from GFF files ===")

    gff37a_path = fp("GCF_000001405.25_GRCh37.p13_genomic.gff")
    gff37b_path = fp("ref_GRCh37.p10_top_level.gff3")
    gff38_path  = fp("GCF_000001405.40_GRCh38.p14_genomic.gff")
    if not os.path.exists(gff38_path): gff38_path += ".gz"

    exact37a, base37a = build_nm_sets(gff37a_path, "GCF_000001405.25 (GRCh37)")
    exact37b, base37b = build_nm_sets(gff37b_path, "ref_GRCh37.p10 (GRCh37)")
    exact38,  base38  = build_nm_sets(gff38_path,  "GCF_000001405.40 (GRCh38)")

    # Combined GRCh37
    exact37 = exact37a | exact37b
    base37  = base37a  | base37b

    print("\n=== Collecting RefSeq IDs to check ===")
    records = collect_ids()
    print(f"  Total IDs collected: {len(records)}")

    print("\n=== Results: IDs NOT found (even base ID absent) ===\n")

    missing_rows = []
    for source, gene, nm, assembly in records:
        nm_base = nm.split(".")[0]
        if assembly == "GRCh37":
            exact_found = nm in exact37
            base_found  = nm_base in base37
        else:
            exact_found = nm in exact38
            base_found  = nm_base in base38

        if not base_found:
            missing_rows.append((source, gene, nm, assembly, "BASE MISSING"))
        elif not exact_found:
            missing_rows.append((source, gene, nm, assembly, "version only"))

    # Separate into two groups for display
    base_missing   = [(s,g,n,a,r) for s,g,n,a,r in missing_rows if r == "BASE MISSING"]
    ver_only       = [(s,g,n,a,r) for s,g,n,a,r in missing_rows if r == "version only"]

    print(f"{'Source':<42} {'Gene':<20} {'NM ID':<20} {'Assembly':<8} {'Status'}")
    print("-" * 110)
    if base_missing:
        print("\n### Base NM ID not found in ANY GFF (truly missing) ###\n")
        for s,g,n,a,r in sorted(base_missing, key=lambda x: (x[3],x[1])):
            print(f"  {s:<40} {g:<20} {n:<20} {a:<8} {r}")
    else:
        print("  (none — all base IDs are present)")

    if ver_only:
        print(f"\n### Specific version not found (but base ID exists in GFF) — {len(ver_only)} entries ###\n")
        for s,g,n,a,r in sorted(ver_only, key=lambda x: (x[3],x[1])):
            print(f"  {s:<40} {g:<20} {n:<20} {a:<8} {r}")
    else:
        print("\n  (none — all versioned IDs found exactly)")

    print(f"\nSummary: {len(base_missing)} truly missing IDs | {len(ver_only)} version-only mismatches | {len(records)-len(missing_rows)} fully found")

if __name__ == "__main__":
    main()
