#!/usr/bin/env python3
"""
Build gene_data.json for the transcript diff visualization web app.

Source of truth for transcript lists:
  files/output/ensembl_transcripts_grch37.tsv
  files/output/ensembl_transcripts_grch38.tsv
  files/output/refseq_transcripts_grch37.tsv
  files/output/refseq_transcripts_grch38.tsv

Source tags are derived from isoform reference files.
Default selection uses priority-based logic, then longest transcript.

Run from the project root:
    python3 pipeline/build_data.py
"""
import re
import os
import sys
import gzip
import json
import difflib
import datetime
from collections import defaultdict
from typing import Dict, Set, Optional, Tuple, List

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES     = os.path.join(ROOT, "files")
OUTPUT    = os.path.join(ROOT, "public", "data", "gene_data.json")

def fp(*parts):       return os.path.join(FILES, *parts)
def finput(*parts):   return os.path.join(FILES, "input", *parts)
def fisoform(*parts): return os.path.join(FILES, "isoform", *parts)
def foutput(*parts):  return os.path.join(FILES, "output", *parts)

# Gene symbol normalisation map (old -> new)
GENE_SYMBOL_UPDATES = {
    "EIF2C1": "AGO1",
    "FTSJD1": "CMTR2",
    "PAK7":   "PAK5",
    "WHSC1":  "NSD2",
    "STK19":  "WHR1",
}

# ---------------------------------------------------------------------------
# Source priority configuration
# ---------------------------------------------------------------------------
# Display order in dropdown (most important first)
SOURCE_DISPLAY_ORDER = [
    "MANE", "MANE_37_lifted", "clinical",
    "mskcc_isoform_37", "mskcc_isoform_38",
    "oncokb_37", "oncokb_38",
]

# Default selection priority per column type
PRIORITY_GRCH37_ENST = ["mskcc_isoform_37", "oncokb_37", "MANE_37_lifted"]
PRIORITY_GRCH37_NM   = ["clinical", "mskcc_isoform_37", "oncokb_37"]
PRIORITY_GRCH38_ENST = ["MANE", "oncokb_38", "mskcc_isoform_38"]
PRIORITY_GRCH38_NM   = ["MANE", "oncokb_38", "mskcc_isoform_38"]

# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------
def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")

# ---------------------------------------------------------------------------
# ID helpers
# ---------------------------------------------------------------------------
def _nm_num(nm: str) -> int:
    m = re.search(r"NM_(\d+)", nm)
    return int(m.group(1)) if m else 999999999

def _enst_num(enst: str) -> int:
    m = re.search(r"ENST(\d+)", enst)
    return int(m.group(1)) if m else 999999999

def _enst_ver(enst: str) -> int:
    m = re.search(r"\.(\d+)$", enst)
    return int(m.group(1)) if m else 0

def _id_ver(tid: str) -> int:
    """Extract version number from any versioned transcript ID."""
    m = re.search(r"\.(\d+)$", tid)
    return int(m.group(1)) if m else 0

def _id_base(tid: str) -> str:
    return tid.split(".")[0] if tid else ""

# ---------------------------------------------------------------------------
# Source formatting
# ---------------------------------------------------------------------------
def sort_sources(sources: Set[str]) -> str:
    """Format source set in priority display order."""
    ordered   = [s for s in SOURCE_DISPLAY_ORDER if s in sources]
    remaining = sorted(s for s in sources if s not in SOURCE_DISPLAY_ORDER)
    return "; ".join(ordered + remaining)

# ---------------------------------------------------------------------------
# Output TSV loader (primary transcript source of truth)
# ---------------------------------------------------------------------------
def load_output_tsv(
    path: str,
    id_col: str,
) -> Tuple[Dict[str, List[Tuple[str, str, int]]], Dict[str, str], Dict[str, str]]:
    """
    Load one of the four output TSVs.

    Returns:
        gene_to_transcripts  {gene_symbol: [(transcript_id, protein_seq, aa_len)]}
        id_to_seq            {transcript_id: protein_seq}
        id_to_gene_id        {transcript_id: ensembl_gene_id}
    """
    gene_to_transcripts: Dict[str, List[Tuple[str, str, int]]] = defaultdict(list)
    id_to_seq:           Dict[str, str] = {}
    id_to_gene_id:       Dict[str, str] = {}

    with open(path, encoding="utf-8") as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        gene_col    = idx["gene_symbol"]
        id_col_i    = idx[id_col]
        seq_col     = idx["protein_sequence"]
        gid_col     = idx.get("ensembl_gene_id", -1)

        seen: Set[Tuple[str, str]] = set()   # (gene, tid) dedup
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(gene_col, id_col_i, seq_col):
                continue
            gene = GENE_SYMBOL_UPDATES.get(fields[gene_col].strip(), fields[gene_col].strip())
            tid  = fields[id_col_i].strip()
            seq  = fields[seq_col].strip()
            gid  = fields[gid_col].strip() if gid_col >= 0 and gid_col < len(fields) else ""

            if not (gene and tid and seq):
                continue
            key = (gene, tid)
            if key in seen:
                continue
            seen.add(key)

            gene_to_transcripts[gene].append((tid, seq, len(seq)))
            id_to_seq[tid]      = seq
            if gid:
                id_to_gene_id[tid] = gid

    return dict(gene_to_transcripts), id_to_seq, id_to_gene_id

# ---------------------------------------------------------------------------
# Source map builder
# ---------------------------------------------------------------------------
def _src_lookup(src_map: Dict[str, Set[str]], tid: str) -> Set[str]:
    """Look up sources for a transcript ID.
    Checks both the full versioned ID and the base ID (for source files that omit versions).
    """
    return src_map.get(tid, set()) | src_map.get(_id_base(tid), set())


def build_source_maps(isoform_dir: str, nm37_ids: Set[str]) -> Dict[str, Dict[str, Dict[str, Set[str]]]]:
    """
    Build source-tag lookup from all isoform reference files.

    Returns {
        'grch37_enst': {gene: {versioned_id_or_base: set(source_names)}},
        'grch37_nm':   {gene: {versioned_id_or_base: set(source_names)}},
        'grch38_enst': {gene: {versioned_id_or_base: set(source_names)}},
        'grch38_nm':   {gene: {versioned_id_or_base: set(source_names)}},
    }
    Keys preserve the version from the source file; use _src_lookup() to query.

    nm37_ids: full set of NM IDs present in refseq_transcripts_grch37.tsv —
              used to guard MANE tagging for GRCh37 (only exact-match IDs are tagged).
    """
    grch37_enst: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    grch37_nm:   Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    grch38_enst: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    grch38_nm:   Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))

    def ng(g: str) -> str:
        return GENE_SYMBOL_UPDATES.get(g.strip(), g.strip())

    # 1. MANE_GRCh37_list_filtered.csv
    #    "Ensembl StableID GRCh37 (Not MANE)" → grch37 ensembl → MANE_37_lifted
    path = os.path.join(isoform_dir, "MANE_GRCh37_list_filtered.csv")
    with _open(path) as fh:
        header_line = fh.readline().rstrip("\n")
        sep  = "\t" if "\t" in header_line else ","
        cols = [c.strip() for c in header_line.split(sep)]
        for line in fh:
            fields = [f.strip() for f in line.rstrip("\n").split(sep)]
            row  = dict(zip(cols, fields))
            gene = ng(row.get("Gene", ""))
            e37  = row.get("Ensembl StableID GRCh37 (Not MANE)", "").strip()
            if gene and e37:
                grch37_enst[gene][e37].add("MANE_37_lifted")

    # 2. oncokb_isoform_versioned.tsv
    #    GRCh37 Isoform → grch37 ensembl → oncokb_37
    #    GRCh37 RefSeq  → grch37 refseq  → oncokb_37
    #    GRCh38 Isoform → grch38 ensembl → oncokb_38
    #    GRCh38 RefSeq  → grch38 refseq  → oncokb_38
    path = os.path.join(isoform_dir, "oncokb_isoform_versioned.tsv")
    with _open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            row  = dict(zip(cols, fields))
            gene = ng(row.get("Hugo Symbol", ""))
            e37  = row.get("GRCh37 Isoform", "").strip()
            n37  = row.get("GRCh37 RefSeq", "").strip()
            e38  = row.get("GRCh38 Isoform", "").strip()
            n38  = row.get("GRCh38 RefSeq", "").strip()
            if not gene:
                continue
            if e37: grch37_enst[gene][e37].add("oncokb_37")
            if n37: grch37_nm[gene][n37].add("oncokb_37")
            if e38: grch38_enst[gene][e38].add("oncokb_38")
            if n38: grch38_nm[gene][n38].add("oncokb_38")

    # 3. MANE.GRCh38.v1.2.summary.txt
    #    RefSeq_nuc  → grch38 + grch37 refseq  → MANE  (NM IDs are assembly-agnostic)
    #    Ensembl_nuc → grch38 ensembl           → MANE
    path = os.path.join(isoform_dir, "MANE.GRCh38.v1.2.summary.txt")
    if not os.path.exists(path):
        path += ".gz"
    with _open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            row  = dict(zip(cols, fields))
            gene = ng(row.get("symbol", ""))
            nm   = row.get("RefSeq_nuc", "").strip()
            enst = row.get("Ensembl_nuc", "").strip()
            if not gene:
                continue
            if nm:
                grch38_nm[gene][nm].add("MANE")
                if nm in nm37_ids:  # only tag GRCh37 if exact versioned ID exists there
                    grch37_nm[gene][nm].add("MANE")
            if enst: grch38_enst[gene][enst].add("MANE")

    # 4. isoform_overrides_at_mskcc_grch37.txt
    #    enst_id    → grch37 ensembl → mskcc_isoform_37
    #    refseq_id  → grch37 refseq  → mskcc_isoform_37
    path = os.path.join(isoform_dir, "isoform_overrides_at_mskcc_grch37.txt")
    with _open(path) as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gene = ng(parts[0])
            nm   = parts[1].strip()
            enst = parts[2].strip()
            if not gene:
                continue
            if enst: grch37_enst[gene][enst].add("mskcc_isoform_37")
            if nm:   grch37_nm[gene][nm].add("mskcc_isoform_37")

    # 5. isoform_overrides_at_mskcc_grch38.txt
    #    enst_id    → grch38 ensembl → mskcc_isoform_38
    #    refseq_id  → grch38 refseq  → mskcc_isoform_38
    path = os.path.join(isoform_dir, "isoform_overrides_at_mskcc_grch38.txt")
    with _open(path) as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gene = ng(parts[0])
            nm   = parts[1].strip()
            enst = parts[2].strip()
            if not gene:
                continue
            if enst: grch38_enst[gene][enst].add("mskcc_isoform_38")
            if nm:   grch38_nm[gene][nm].add("mskcc_isoform_38")

    # 6. Iv7_dmp_isoform_merged_overrides.txt
    #    "#isoform_override" where value starts with NM_ → grch37 refseq → clinical
    path = os.path.join(isoform_dir, "Iv7_dmp_isoform_merged_overrides.txt")
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            isoform_override = parts[0].strip()
            gene = ng(parts[1])
            if gene and isoform_override.startswith("NM_"):
                if isoform_override:
                    grch37_nm[gene][isoform_override].add("clinical")

    return {
        "grch37_enst": {g: dict(v) for g, v in grch37_enst.items()},
        "grch37_nm":   {g: dict(v) for g, v in grch37_nm.items()},
        "grch38_enst": {g: dict(v) for g, v in grch38_enst.items()},
        "grch38_nm":   {g: dict(v) for g, v in grch38_nm.items()},
    }

# ---------------------------------------------------------------------------
# Default transcript selection
# ---------------------------------------------------------------------------
def select_default(
    candidates: List[Tuple[str, str, int]],   # (tid, seq, aa_len)
    source_map: Dict[str, Set[str]],           # id_base -> set(sources)
    priority: List[str],
) -> str:
    """
    Pick the default transcript:
      1. Highest-priority source group (from priority list)
      2. Within that group, longest protein sequence
      3. If no priority source matches, pick the longest overall
    """
    if not candidates:
        return ""
    # Sort key: (aa_length DESC, version DESC) — longer and newer wins on ties
    sort_key = lambda x: (x[2], _id_ver(x[0]))
    for src in priority:
        matches = [
            (tid, seq, ln) for tid, seq, ln in candidates
            if src in _src_lookup(source_map, tid)
        ]
        if matches:
            return max(matches, key=sort_key)[0]
    # No priority match — longest overall, then highest version
    return max(candidates, key=sort_key)[0]

# ---------------------------------------------------------------------------
# MANE parsing
# ---------------------------------------------------------------------------
def parse_mane_grch38(path: str) -> Dict[str, Dict]:
    """Returns {ENST.v -> {nm, status}}"""
    result: Dict[str, Dict] = {}
    with _open(path) as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = {h.lower(): i for i, h in enumerate(header)}
        enst_col   = next((idx[k] for k in idx if "ensembl" in k and "nuc" in k), None)
        nm_col     = next((idx[k] for k in idx if "refseq"  in k and "nuc" in k), None)
        status_col = next((idx[k] for k in idx if "mane_status" in k or k == "status"), None)
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if enst_col is None or enst_col >= len(fields):
                continue
            enst   = fields[enst_col].strip()
            nm     = fields[nm_col].strip()     if nm_col     is not None and nm_col     < len(fields) else ""
            status = fields[status_col].strip() if status_col is not None and status_col < len(fields) else ""
            if enst.startswith("ENST"):
                result[enst] = {"nm": nm, "status": status}
    return result

def parse_mane_grch37_csv(path: str) -> Dict[str, Dict]:
    """
    Returns {gene -> {grch38_enst, grch37_enst, nm, mane_type}}
    """
    result: Dict[str, Dict] = {}
    with _open(path) as fh:
        header_line = fh.readline().rstrip("\n")
        sep  = "\t" if "\t" in header_line else ","
        cols = header_line.split(sep)
        for line in fh:
            fields = line.rstrip("\n").split(sep)
            if len(fields) < 5:
                continue
            row = dict(zip(cols, fields))
            gene        = GENE_SYMBOL_UPDATES.get(row.get("Gene", "").strip(), row.get("Gene", "").strip())
            mane_type   = row.get("MANE TYPE", "").strip()
            grch38_enst = row.get("Ensembl StableID GRCh38", "").strip()
            nm          = row.get("RefSeq StableID GRCh38 / GRCh37", "").strip()
            grch37_enst = row.get("Ensembl StableID GRCh37 (Not MANE)", "").strip()
            if gene:
                result[gene] = {
                    "grch38_enst": grch38_enst,
                    "grch37_enst": grch37_enst,
                    "nm":          nm,
                    "mane_type":   mane_type,
                }
    return result

# ---------------------------------------------------------------------------
# Clinical / annotation file parsing
# ---------------------------------------------------------------------------
def parse_iv7_overrides(path: str) -> Dict[str, str]:
    """Returns {gene -> clinical_refseq base NM}"""
    result: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            gene = GENE_SYMBOL_UPDATES.get(parts[1].strip(), parts[1].strip())
            nm   = parts[2].strip()
            if gene and nm.startswith("NM_"):
                result[gene] = nm.split(".")[0]
    return result

def parse_iv7_isoform_override(path: str) -> Dict[str, str]:
    """Returns {gene -> isoform_override}"""
    result: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            isoform_override = parts[0].strip()
            gene = GENE_SYMBOL_UPDATES.get(parts[1].strip(), parts[1].strip())
            if gene and isoform_override:
                result[gene] = isoform_override
    return result

def parse_mskcc_overrides(path: str) -> Dict[str, Dict]:
    """
    isoform_overrides_at_mskcc_grch37.txt
    Returns {gene -> {refseq_id, enst_id, note}}
    """
    result: Dict[str, Dict] = {}
    with open(path, "r", encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gene   = GENE_SYMBOL_UPDATES.get(parts[0].strip(), parts[0].strip())
            refseq = parts[1].strip()
            enst   = parts[2].strip()
            note   = parts[3].strip() if len(parts) > 3 else ""
            if gene:
                result[gene] = {"refseq_id": refseq, "enst_id": enst, "note": note}
    return result

def parse_mskcc_grch38_overrides(path: str) -> Dict[str, Dict]:
    """
    isoform_overrides_at_mskcc_grch38.txt
    Returns {gene -> {refseq_id, enst_id, note}}
    """
    result: Dict[str, Dict] = {}
    with open(path, "r", encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gene   = GENE_SYMBOL_UPDATES.get(parts[0].strip(), parts[0].strip())
            refseq = parts[1].strip()
            enst   = parts[2].strip()
            note   = parts[3].strip() if len(parts) > 3 else ""
            if gene:
                result[gene] = {"refseq_id": refseq, "enst_id": enst, "note": note}
    return result

def parse_germline(path: str) -> Set[str]:
    genes: Set[str] = set()
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            g = line.strip()
            if g:
                genes.add(GENE_SYMBOL_UPDATES.get(g, g))
    return genes

def parse_oncokb_isoform(path: str) -> Dict[str, Dict]:
    """
    oncokb_isoform_versioned.tsv
    Returns {gene -> {entrez_gene_id, grch37_enst_base, grch37_nm, grch38_enst_base, grch38_nm, gene_type}}
    """
    result: Dict[str, Dict] = {}
    with open(path, "r", encoding="utf-8") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            row  = dict(zip(cols, fields))
            gene = row.get("Hugo Symbol", "").strip()
            gene = re.sub(r'\s+\(', '(', gene)
            gene = GENE_SYMBOL_UPDATES.get(gene, gene)
            if not gene:
                continue
            result[gene] = {
                "entrez_gene_id":  row.get("Entrez Gene ID", "").strip(),
                "grch37_enst_base": row.get("GRCh37 Isoform", "").strip(),
                "grch37_nm":       row.get("GRCh37 RefSeq", "").strip(),
                "grch38_enst_base": row.get("GRCh38 Isoform", "").strip(),
                "grch38_nm":       row.get("GRCh38 RefSeq", "").strip(),
                "gene_type":       row.get("Gene Type", "").strip(),
            }
    return result

def parse_hgnc(path: str) -> Tuple[Dict[str, Dict], Dict[str, str], Dict[str, str]]:
    """
    Parse HGNC complete set TSV.
    Returns (hgnc_map, alias_map, prev_map)
    """
    hgnc_map:  Dict[str, Dict] = {}
    alias_map: Dict[str, str]  = {}
    prev_map:  Dict[str, str]  = {}
    with _open(path) as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = {h.strip(): i for i, h in enumerate(header)}
        hgnc_col   = idx.get("hgnc_id", 0)
        symbol_col = idx.get("symbol", 1)
        alias_col  = idx.get("alias_symbol", -1)
        prev_col   = idx.get("prev_symbol", 10)
        entrez_col = idx.get("entrez_id", 18)
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(hgnc_col, symbol_col, entrez_col):
                continue
            hgnc_id = fields[hgnc_col].strip()
            symbol  = fields[symbol_col].strip()
            raw_eid = fields[entrez_col].strip()
            entrez_id = ""
            if raw_eid and raw_eid != "nan":
                try:
                    entrez_id = str(int(float(raw_eid)))
                except (ValueError, OverflowError):
                    entrez_id = raw_eid
            info = {"hgnc_id": hgnc_id, "entrez_id": entrez_id}
            if symbol:
                hgnc_map[symbol] = info
            if alias_col >= 0 and alias_col < len(fields):
                raw_alias = fields[alias_col].strip()
                if raw_alias and raw_alias != "nan":
                    for a in raw_alias.split("|"):
                        a = a.strip().strip('"')
                        if a and a not in hgnc_map and a not in alias_map:
                            alias_map[a] = symbol
            if prev_col < len(fields):
                raw_prev = fields[prev_col].strip()
                if raw_prev and raw_prev != "nan":
                    for ps in raw_prev.split("|"):
                        ps = ps.strip().strip('"')
                        if ps and ps not in hgnc_map and ps not in prev_map:
                            prev_map[ps] = symbol
                            hgnc_map.setdefault(ps, info)
    return hgnc_map, alias_map, prev_map

# ---------------------------------------------------------------------------
# Similarity
# ---------------------------------------------------------------------------
def compute_similarity(seq1: str, seq2: str) -> Dict:
    if not seq1 or not seq2:
        return {"pct": None, "diff_count": 0}
    sm = difflib.SequenceMatcher(None, seq1, seq2, autojunk=False)
    pct = round(sm.ratio() * 100, 4)
    diff_count = sum(
        max(i2 - i1, j2 - j1)
        for tag, i1, i2, j1, j2 in sm.get_opcodes()
        if tag != "equal"
    )
    return {"pct": pct, "diff_count": diff_count}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=== Transcript Diff Pipeline ===")

    print("\n[1/8] Loading transcript output TSVs (source of truth)…")
    gene_to_enst37, enst37_to_seq, enst37_to_gene_id = load_output_tsv(
        foutput("ensembl_transcripts_grch37.tsv"), "ensembl_transcript_id"
    )
    gene_to_enst38, enst38_to_seq, enst38_to_gene_id = load_output_tsv(
        foutput("ensembl_transcripts_grch38.tsv"), "ensembl_transcript_id"
    )
    gene_to_nm37, nm37_to_seq, _ = load_output_tsv(
        foutput("refseq_transcripts_grch37.tsv"), "refseq_transcript_id"
    )
    gene_to_nm38, nm38_to_seq, _ = load_output_tsv(
        foutput("refseq_transcripts_grch38.tsv"), "refseq_transcript_id"
    )
    # Build reverse: versioned ENST -> gene symbol (from GRCh38)
    enst38_to_gene_sym: Dict[str, str] = {}
    for g, txs in gene_to_enst38.items():
        for tid, seq, ln in txs:
            enst38_to_gene_sym[tid] = g

    n37_genes = sum(len(v) for v in gene_to_enst37.values())
    n38_genes = sum(len(v) for v in gene_to_enst38.values())
    print(f"  GRCh37: {n37_genes:,} ENSTs, {sum(len(v) for v in gene_to_nm37.values()):,} NMs")
    print(f"  GRCh38: {n38_genes:,} ENSTs, {sum(len(v) for v in gene_to_nm38.values()):,} NMs")

    print("[2/8] Parsing MANE files…")
    mane38_path = fisoform("MANE.GRCh38.v1.2.summary.txt")
    if not os.path.exists(mane38_path):
        mane38_path += ".gz"
    mane_grch38     = parse_mane_grch38(mane38_path)   # {ENST.v -> {nm, status}}
    mane_grch37_csv = parse_mane_grch37_csv(fisoform("MANE_GRCh37_list_filtered.csv"))
    print(f"  {len(mane_grch38):,} MANE GRCh38, {len(mane_grch37_csv):,} MANE GRCh37 entries")

    print("[3/8] Parsing clinical / annotation files…")
    iv7          = parse_iv7_overrides(fisoform("Iv7_dmp_isoform_merged_overrides.txt"))
    iv7_override = parse_iv7_isoform_override(fisoform("Iv7_dmp_isoform_merged_overrides.txt"))
    mskcc        = parse_mskcc_overrides(fisoform("isoform_overrides_at_mskcc_grch37.txt"))
    mskcc38      = parse_mskcc_grch38_overrides(fisoform("isoform_overrides_at_mskcc_grch38.txt"))
    germline     = parse_germline(fisoform("germline_panel_94.txt"))
    print(f"  Iv7: {len(iv7)}, MSKCC37: {len(mskcc)}, MSKCC38: {len(mskcc38)}, germline: {len(germline)}")

    print("[4/8] Parsing OncoKB isoform file…")
    oncokb_isoform = parse_oncokb_isoform(fisoform("oncokb_isoform_versioned.tsv"))
    print(f"  OncoKB isoform entries: {len(oncokb_isoform)}")

    print("[5/8] Parsing HGNC complete set…")
    hgnc_map, alias_map, prev_map = parse_hgnc(fisoform("hgnc_complete_set_oct_07_2025.txt"))
    print(f"  HGNC entries: {len(hgnc_map)}")

    print("[6/8] Building source maps from isoform reference files…")
    src_maps = build_source_maps(fisoform(), nm37_ids=set(nm37_to_seq.keys()))
    total_tags = sum(
        len(v) for d in src_maps.values() for v in d.values()
    )
    print(f"  Source tags built: {total_tags:,} (gene, id, source) combinations")

    all_genes: Set[str] = set(iv7) | set(mskcc) | set(mskcc38) | set(oncokb_isoform)
    print(f"\n[7/8] Building gene records for {len(all_genes)} genes…")

    genes_out: List[Dict] = []
    warnings:  List[str]  = []

    for gene in sorted(all_genes):
        mskcc_entry   = mskcc.get(gene, {})
        mskcc38_entry = mskcc38.get(gene, {})
        iv7_nm_base   = iv7.get(gene)
        is_germline   = gene in germline
        is_clinical   = gene in iv7
        is_oncokb     = gene in oncokb_isoform
        okb_iso       = oncokb_isoform.get(gene, {})

        hgnc_entry     = hgnc_map.get(gene, {})
        hgnc_id        = hgnc_entry.get("hgnc_id", "")
        entrez_gene_id = hgnc_entry.get("entrez_id", "") or okb_iso.get("entrez_gene_id", "")

        # ---- Transcript lists from output TSVs ----
        e37 = gene_to_enst37.get(gene, [])   # [(tid, seq, aa_len)]
        e38 = gene_to_enst38.get(gene, [])
        n37 = gene_to_nm37.get(gene, [])
        n38 = gene_to_nm38.get(gene, [])

        # ---- Source maps for this gene ----
        e37_src = src_maps["grch37_enst"].get(gene, {})   # {id_base: set(sources)}
        e38_src = src_maps["grch38_enst"].get(gene, {})
        n37_src = src_maps["grch37_nm"].get(gene, {})
        n38_src = src_maps["grch38_nm"].get(gene, {})

        # ---- Default selection: priority source, then longest ----
        grch37_enst = select_default(e37, e37_src, PRIORITY_GRCH37_ENST)
        grch38_enst = select_default(e38, e38_src, PRIORITY_GRCH38_ENST)
        grch37_nm   = select_default(n37, n37_src, PRIORITY_GRCH37_NM)
        grch38_nm   = select_default(n38, n38_src, PRIORITY_GRCH38_NM)

        # ---- Build transcript objects ----
        all_transcripts: List[Dict] = []

        for tid, seq, ln in e37:
            sources = _src_lookup(e37_src, tid)
            all_transcripts.append({
                "id":         tid,
                "type":       "ensembl",
                "assembly":   "GRCh37",
                "sequence":   seq,
                "length":     ln,
                "source":     sort_sources(sources),
                "is_primary": tid == grch37_enst,
            })
        for tid, seq, ln in e38:
            sources = _src_lookup(e38_src, tid)
            all_transcripts.append({
                "id":         tid,
                "type":       "ensembl",
                "assembly":   "GRCh38",
                "sequence":   seq,
                "length":     ln,
                "source":     sort_sources(sources),
                "is_primary": tid == grch38_enst,
            })
        for tid, seq, ln in n37:
            sources = _src_lookup(n37_src, tid)
            all_transcripts.append({
                "id":         tid,
                "type":       "refseq",
                "assembly":   "GRCh37",
                "sequence":   seq,
                "length":     ln,
                "source":     sort_sources(sources),
                "is_primary": tid == grch37_nm,
            })
        for tid, seq, ln in n38:
            sources = _src_lookup(n38_src, tid)
            all_transcripts.append({
                "id":         tid,
                "type":       "refseq",
                "assembly":   "GRCh38",
                "sequence":   seq,
                "length":     ln,
                "source":     sort_sources(sources),
                "is_primary": tid == grch38_nm,
            })

        # ---- MANE GRCh38 info ----
        mane38_enst = ""
        mane38_nm   = ""
        mane38_status = ""
        # Check selected grch38_enst first
        if grch38_enst and grch38_enst in mane_grch38:
            mane38_enst   = grch38_enst
            mane38_nm     = mane_grch38[grch38_enst].get("nm", "")
            mane38_status = mane_grch38[grch38_enst].get("status", "")
        # If not, scan gene's GRCh38 transcripts for MANE Select
        if not mane38_enst:
            for tid, seq, ln in e38:
                s = mane_grch38.get(tid, {}).get("status", "")
                if "MANE Select" in s:
                    mane38_enst   = tid
                    mane38_nm     = mane_grch38[tid].get("nm", "")
                    mane38_status = s
                    break

        # ---- MANE GRCh37 validated ----
        mane37_entry       = mane_grch37_csv.get(gene, {})
        mane37_grch38_enst = mane37_entry.get("grch38_enst", "")
        mane37_grch37_enst = mane37_entry.get("grch37_enst", "")
        mane_grch37_valid  = mane37_grch38_enst in mane_grch38 if mane37_grch38_enst else False

        curated_note = mskcc_entry.get("note", "")

        # MSKCC GRCh38 IDs (as-is; UI can display them in the collections)
        mskcc38_enst_raw = mskcc38_entry.get("enst_id", "")
        mskcc38_nm_raw   = mskcc38_entry.get("refseq_id", "")

        # mane_only flag: True when gene has no curated entry except MANE
        _has_curated = bool(
            mskcc_entry or mskcc38_entry or iv7_nm_base or
            okb_iso.get("grch37_enst_base") or okb_iso.get("grch38_enst_base")
        )
        mane_only = not _has_curated

        # ---- Sequences for similarity ----
        seq_map = {t["id"]: t["sequence"] for t in all_transcripts}
        s37e   = seq_map.get(grch37_enst, "")
        s38e   = seq_map.get(grch38_enst, "")
        s37n   = seq_map.get(grch37_nm, "")
        s38n   = seq_map.get(grch38_nm, "")
        s_mane = seq_map.get(mane38_enst, "")

        similarities = {
            "grch37_enst_vs_grch37_refseq": compute_similarity(s37e, s37n),
            "grch38_enst_vs_grch38_refseq": compute_similarity(s38e, s38n),
            "grch37_enst_vs_grch38_enst":   compute_similarity(s37e, s38e),
            "grch37_enst_vs_mane":          compute_similarity(s37e, s_mane),
            "grch38_enst_vs_mane":          compute_similarity(s38e, s_mane),
            "grch37_nm_vs_mane":            compute_similarity(s37n, s_mane),
        }

        # ---- Ensembl Gene ID (from TSV) ----
        gene_id = ""
        if grch37_enst:
            gene_id = enst37_to_gene_id.get(grch37_enst, "")
        if not gene_id and grch38_enst:
            gene_id = enst38_to_gene_id.get(grch38_enst, "")

        # ---- Warnings ----
        if is_clinical and not grch37_nm:
            warnings.append(f"NO_GRCh37_NM {gene}: clinical gene has no GRCh37 NM in output TSV")
        if is_clinical and not grch38_nm:
            warnings.append(f"NO_GRCh38_NM {gene}: clinical gene has no GRCh38 NM in output TSV")
        if is_clinical and not grch37_enst:
            warnings.append(f"NO_GRCh37_ENST {gene}: clinical gene has no GRCh37 ENST in output TSV")
        if is_germline and s37e and s38e and s37e != s38e:
            warnings.append(f"GERMLINE_SEQ_MISMATCH {gene}: GRCh37 vs GRCh38 ENST sequences differ")

        gene_record = {
            "gene_symbol":     gene,
            "gene_id":         gene_id,
            "hgnc_id":         hgnc_id,
            "entrez_gene_id":  entrez_gene_id,
            "is_germline":     is_germline,
            "is_clinical":     is_clinical,
            "is_oncokb":       is_oncokb,
            "mane_only":       mane_only,
            "clinical_refseq": iv7_nm_base or "",
            "curated_note":    curated_note,
            "mskcc38_enst":    mskcc38_enst_raw or "",
            "mskcc38_nm":      mskcc38_nm_raw or "",
            "mane_grch38": {
                "enst":   mane38_enst or "",
                "nm":     mane38_nm or "",
                "status": mane38_status or "",
            },
            "mane_grch37_enst": mane37_grch37_enst if mane_grch37_valid else "",
            "best_match": {
                "grch37_enst": grch37_enst or "",
                "grch38_enst": grch38_enst or "",
                "grch37_nm":   grch37_nm or "",
                "grch38_nm":   grch38_nm or "",
            },
            "transcripts":  all_transcripts,
            "similarities": similarities,
        }
        genes_out.append(gene_record)

    # ---- Collections ----
    print("[8/8] Building collections…")
    collections: Dict = {
        "oncokb": {
            g: {"grch38_enst": v.get("grch38_enst_base", ""), "grch38_nm": v.get("grch38_nm", "")}
            for g, v in oncokb_isoform.items()
        },
        "oncokb_isoform": {
            g: {
                "entrez_gene_id": v.get("entrez_gene_id", ""),
                "grch37_enst":    v.get("grch37_enst_base", ""),
                "grch37_nm":      v.get("grch37_nm", ""),
                "grch38_enst":    v.get("grch38_enst_base", ""),
                "grch38_nm":      v.get("grch38_nm", ""),
                "gene_type":      v.get("gene_type", ""),
            }
            for g, v in oncokb_isoform.items()
        },
        "clinical": {
            g: {"grch37_nm": nm, "isoform_override": iv7_override.get(g, "")}
            for g, nm in iv7.items()
        },
        "mane": {},
        "mskcc_isoform": {
            g: {"grch37_enst": v.get("enst_id", ""), "grch37_nm": v.get("refseq_id", "")}
            for g, v in mskcc.items()
        },
        "mskcc_grch38_isoform": {
            g: {"grch38_enst": v.get("enst_id", ""), "grch38_nm": v.get("refseq_id", "")}
            for g, v in mskcc38.items()
        },
    }
    # Build mane collection from mane_grch38 (ENST.v -> {nm, status})
    for enst, info in mane_grch38.items():
        gene_sym = enst38_to_gene_sym.get(enst, "")
        if gene_sym:
            collections["mane"][gene_sym] = {
                "grch38_enst": enst,
                "grch38_nm":   info.get("nm", ""),
                "status":      info.get("status", ""),
            }

    output = {
        "metadata": {
            "generated_at":    datetime.datetime.utcnow().isoformat() + "Z",
            "ensembl_version": 111,
            "mane_version":    "1.2",
            "gene_count":      len(genes_out),
        },
        "genes":       genes_out,
        "collections": collections,
    }

    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    with open(OUTPUT, "w", encoding="utf-8") as fh:
        json.dump(output, fh, separators=(",", ":"))

    size_mb = os.path.getsize(OUTPUT) / 1_048_576
    print(f"\nWrote {len(genes_out)} genes → {OUTPUT} ({size_mb:.1f} MB)")

    if warnings:
        print(f"\n=== {len(warnings)} WARNINGS ===")
        for w in sorted(set(warnings))[:60]:
            print(f"  {w}")
        if len(warnings) > 60:
            print(f"  … and {len(warnings) - 60} more")

if __name__ == "__main__":
    main()
