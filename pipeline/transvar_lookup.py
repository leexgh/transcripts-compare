#!/usr/bin/env python3
"""
TransVar multi-version variant coordinate lookup wrapper.

Usage:
    python3 pipeline/transvar_lookup.py -i "PIK3CA:c.1633G>A"
    python3 pipeline/transvar_lookup.py -l variants.txt --json > results.json

Prerequisites:
    1. Install TransVar:  cd /path/to/transvar && pip install -e .
    2. Download refs:     transvar config --download_ref --refversion hg19
                          transvar config --download_ref --refversion hg38
    3. Index GTF files:   transvar index --ensembl /path/to/GTF.gtf.gz

Output columns (TSV): input | gtf_version | transcript | gene | strand | gDNA | cDNA | protein | region | info
"""
import argparse
import json
import subprocess
import sys
import re
from typing import List, Dict, Optional

# ---------------------------------------------------------------------------
# Configuration: GTF versions and their indexed database paths.
# Edit this list to add/remove GTF versions. Each entry needs:
#   label:     display name for the version
#   refversion: hg19 or hg38
#   db:        path to the TransVar-indexed database (the .transvardb file)
#
# After indexing a GTF with `transvar index --ensembl <gtf>`,
# the database is created at <gtf>.transvardb
# ---------------------------------------------------------------------------
GTF_VERSIONS: List[Dict[str, str]] = [
    # Example entries — uncomment and edit paths after indexing:
    # {
    #     "label": "Ensembl GRCh37.75",
    #     "refversion": "hg19",
    #     "db": "/path/to/files/Homo_sapiens.GRCh37.75.gtf.gz.transvardb",
    # },
    # {
    #     "label": "Ensembl GRCh37.87",
    #     "refversion": "hg19",
    #     "db": "/path/to/files/Homo_sapiens.GRCh37.87.gtf.gz.transvardb",
    # },
    # {
    #     "label": "Ensembl GRCh38.111",
    #     "refversion": "hg38",
    #     "db": "/path/to/files/Homo_sapiens.GRCh38.111.gtf.gz.transvardb",
    # },
]


def detect_anno_type(variant: str) -> str:
    """Auto-detect TransVar annotation type from the variant string."""
    if ":c." in variant:
        return "canno"
    elif ":p." in variant:
        return "panno"
    elif ":g." in variant:
        return "ganno"
    # Default to canno
    return "canno"


def parse_coordinates(coord_field: str) -> Dict[str, str]:
    """
    Parse TransVar coordinates field like 'chr3:g.178936091G>A/c.1633G>A/p.E545K'
    into separate gDNA, cDNA, protein fields.
    """
    parts = coord_field.split("/")
    result = {"gDNA": "", "cDNA": "", "protein": ""}
    for part in parts:
        part = part.strip()
        if ":g." in part or part.startswith("chr"):
            result["gDNA"] = part
        elif ":c." in part or part.startswith("c."):
            result["cDNA"] = part
        elif ":p." in part or part.startswith("p."):
            result["protein"] = part
    return result


def run_transvar(variant: str, anno_type: str, refversion: str, db_path: str) -> List[Dict]:
    """Run TransVar for a single variant against a specific database."""
    cmd = [
        "transvar", anno_type,
        "-i", variant,
        "--refversion", refversion,
        "--ensembl", db_path,
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except FileNotFoundError:
        print("Error: transvar not found. Install it first.", file=sys.stderr)
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print(f"Warning: TransVar timed out for {variant}", file=sys.stderr)
        return []

    results = []
    for line in proc.stdout.strip().split("\n"):
        if not line or line.startswith("input"):
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        # TransVar output: input, transcript, gene, strand, coordinates, region, info
        coords = parse_coordinates(fields[4]) if len(fields) > 4 else {}
        results.append({
            "input": fields[0],
            "transcript": fields[1],
            "gene": fields[2],
            "strand": fields[3],
            "gDNA": coords.get("gDNA", ""),
            "cDNA": coords.get("cDNA", ""),
            "protein": coords.get("protein", ""),
            "region": fields[5] if len(fields) > 5 else "",
            "info": fields[6] if len(fields) > 6 else "",
        })
    return results


def lookup_variant(variant: str, gtf_versions: List[Dict[str, str]]) -> List[Dict]:
    """Run TransVar for a variant across all configured GTF versions."""
    anno_type = detect_anno_type(variant)
    all_results = []
    for ver in gtf_versions:
        results = run_transvar(variant, anno_type, ver["refversion"], ver["db"])
        for r in results:
            r["gtf_version"] = ver["label"]
        all_results.extend(results)
    return all_results


def format_tsv(results: List[Dict]) -> str:
    """Format results as TSV."""
    cols = ["input", "gtf_version", "transcript", "gene", "strand", "gDNA", "cDNA", "protein", "region", "info"]
    lines = ["\t".join(cols)]
    for r in results:
        lines.append("\t".join(r.get(c, "") for c in cols))
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="TransVar multi-version variant lookup")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", help="Single variant (e.g. 'PIK3CA:c.1633G>A')")
    group.add_argument("-l", "--list", help="File with one variant per line")
    parser.add_argument("--json", action="store_true", help="Output JSON instead of TSV")
    args = parser.parse_args()

    if not GTF_VERSIONS:
        print("Error: No GTF versions configured.", file=sys.stderr)
        print("Edit GTF_VERSIONS in pipeline/transvar_lookup.py to add indexed GTF databases.", file=sys.stderr)
        sys.exit(1)

    variants = []
    if args.input:
        variants = [args.input]
    else:
        with open(args.list, "r") as fh:
            variants = [line.strip() for line in fh if line.strip()]

    all_results = []
    for v in variants:
        all_results.extend(lookup_variant(v, GTF_VERSIONS))

    if args.json:
        print(json.dumps(all_results, indent=2))
    else:
        print(format_tsv(all_results))


if __name__ == "__main__":
    main()
