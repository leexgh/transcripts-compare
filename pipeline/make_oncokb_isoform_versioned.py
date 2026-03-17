#!/usr/bin/env python3
"""
Produce oncokb_isoform_versioned.tsv from oncokb_isoform.tsv:

  Step 1: Version ENST IDs
    - GRCh37 Isoform  → look up version in Homo_sapiens.GRCh37.87.gtf
    - GRCh38 Isoform  → look up version in Homo_sapiens.GRCh38.111.gtf

  Step 2: Overwrite GRCh38 Isoform + GRCh38 RefSeq
    - For each gene, if present in oncokb.csv use its ensembl_transcript_id / reference_sequence_id
    - Report genes in oncokb.csv that are NOT in oncokb_isoform.tsv

Run from project root:
    python3 pipeline/make_oncokb_isoform_versioned.py
"""
import re
import os
import gzip

ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(ROOT, "files")

def fp(*parts):
    return os.path.join(FILES, *parts)

def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")

# ---------------------------------------------------------------------------
# 1. Build base_id -> versioned_id maps from GTF files
# ---------------------------------------------------------------------------
def build_enst_version_map(gtf_path: str) -> dict:
    """Returns {ENST_base -> ENST_base.version}"""
    m = {}
    with _open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#") or "\ttranscript\t" not in line:
                continue
            attrs = line.split("\t", 9)[-1] if len(line.split("\t")) >= 9 else ""
            m_id  = re.search(r'transcript_id "([^"]+)"', attrs)
            m_ver = re.search(r'transcript_version "([^"]+)"', attrs)
            if not m_id:
                continue
            base = m_id.group(1)
            ver  = m_ver.group(1) if m_ver else None
            versioned = f"{base}.{ver}" if ver else base
            m[base] = versioned
    return m

def version_enst(enst_base: str, ver_map: dict) -> str:
    """Return versioned ID if found, else original."""
    base = enst_base.split(".")[0]  # strip any existing version
    return ver_map.get(base, enst_base)

# ---------------------------------------------------------------------------
# 2. Parse oncokb.csv  →  {hugo_symbol -> {enst, nm}}  (GRCh38 only)
# ---------------------------------------------------------------------------
def parse_oncokb_csv(path: str) -> dict:
    result = {}
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split(",")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split(",")))
            gene   = row.get("hugo_symbol", "").strip().strip('"')
            genome = row.get("reference_genome", "").strip()
            enst   = row.get("ensembl_transcript_id", "").strip().strip('"')
            nm     = row.get("reference_sequence_id", "").strip().strip('"')
            if gene and genome == "GRCh38":
                result[gene] = {"enst": enst, "nm": nm}
    return result

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Building GRCh37 version map from GTF…")
    gtf37 = fp("Homo_sapiens.GRCh37.87.gtf")
    if not os.path.exists(gtf37): gtf37 += ".gz"
    ver37 = build_enst_version_map(gtf37)
    print(f"  {len(ver37):,} GRCh37 ENST bases")

    print("Building GRCh38 version map from GTF…")
    gtf38 = fp("Homo_sapiens.GRCh38.111.gtf")
    if not os.path.exists(gtf38): gtf38 += ".gz"
    ver38 = build_enst_version_map(gtf38)
    print(f"  {len(ver38):,} GRCh38 ENST bases")

    print("Parsing oncokb.csv…")
    oncokb = parse_oncokb_csv(fp("oncokb.csv"))
    print(f"  {len(oncokb):,} GRCh38 genes in oncokb.csv")

    # Read TSV
    tsv_in  = fp("oncokb_isoform.tsv")
    tsv_out = fp("oncokb_isoform_versioned.tsv")

    with open(tsv_in, encoding="utf-8") as fh:
        lines = fh.readlines()

    header_line = lines[0].rstrip("\n")
    cols = header_line.split("\t")

    # Locate column indices
    def ci(name):
        for i, c in enumerate(cols):
            if c.strip() == name:
                return i
        raise ValueError(f"Column not found: {name!r}")

    idx_gene   = ci("Hugo Symbol")
    idx_37e    = ci("GRCh37 Isoform")
    idx_37r    = ci("GRCh37 RefSeq")
    idx_38e    = ci("GRCh38 Isoform")
    idx_38r    = ci("GRCh38 RefSeq")

    tsv_genes = set()
    out_rows  = [header_line + "\n"]

    for raw in lines[1:]:
        row = raw.rstrip("\n").split("\t")
        if not row or not row[0].strip():
            out_rows.append(raw)
            continue

        gene = row[idx_gene].strip()
        tsv_genes.add(gene)

        # Step 1: version the ENST IDs
        if idx_37e < len(row) and row[idx_37e].startswith("ENST"):
            row[idx_37e] = version_enst(row[idx_37e], ver37)
        if idx_38e < len(row) and row[idx_38e].startswith("ENST"):
            row[idx_38e] = version_enst(row[idx_38e], ver38)

        # Step 2: overwrite GRCh38 columns from oncokb.csv
        if gene in oncokb:
            okb = oncokb[gene]
            if okb["enst"]:
                row[idx_38e] = okb["enst"]
            if okb["nm"]:
                row[idx_38r] = okb["nm"]

        out_rows.append("\t".join(row) + "\n")

    # Report genes in oncokb.csv but not in oncokb_isoform.tsv
    missing = sorted(g for g in oncokb if g not in tsv_genes)
    if missing:
        print(f"\nGenes in oncokb.csv (GRCh38) but NOT in oncokb_isoform.tsv ({len(missing)}):")
        for g in missing:
            okb = oncokb[g]
            print(f"  {g:30s}  ENST={okb['enst']}  NM={okb['nm']}")
    else:
        print("\nAll oncokb.csv genes are present in oncokb_isoform.tsv.")

    with open(tsv_out, "w", encoding="utf-8") as fh:
        fh.writelines(out_rows)

    print(f"\nWrote {len(out_rows)-1} data rows → {tsv_out}")

if __name__ == "__main__":
    main()
