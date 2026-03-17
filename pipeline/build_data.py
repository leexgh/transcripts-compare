#!/usr/bin/env python3
"""
Build gene_data.json for the transcript diff visualization web app.

Run from the project root:
    python3 pipeline/build_data.py
"""
import re
import os
import sys
import gzip
import json
import pickle
import difflib
import datetime
from collections import defaultdict
from typing import Dict, Set, Optional, Tuple, List

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(ROOT, "files")
OUTPUT = os.path.join(ROOT, "public", "data", "gene_data.json")
CACHE_DIR = os.path.join(ROOT, "pipeline", ".cache")
os.makedirs(CACHE_DIR, exist_ok=True)

def fp(*parts):
    return os.path.join(FILES, *parts)

def finput(*parts):
    return os.path.join(FILES, "input", *parts)

def fisoform(*parts):
    return os.path.join(FILES, "isoform", *parts)

# Gene symbol normalisation map (old -> new)
GENE_SYMBOL_UPDATES = {
    "EIF2C1": "AGO1",
    "FTSJD1": "CMTR2",
    "PAK7":   "PAK5",
    "WHSC1":  "NSD2",
    "STK19":  "WHR1",
}

# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------
def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")

def _iter_fasta(path: str):
    """Yield (header, sequence) stripping trailing '*'."""
    with _open(path) as fh:
        header = None
        chunks: List[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).replace("*", "").upper()
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks).replace("*", "").upper()

_RE_ENST = re.compile(r"transcript:(ENST\d+\.\d+)")
_RE_NP   = re.compile(r"^(NP_\d+\.\d+)")

def parse_ensembl_pep(path: str):
    """Returns seq_to_ensts: {seq -> set(ENST.v)}, enst_to_seq: {ENST.v -> seq}"""
    seq_to_ensts: Dict[str, Set[str]] = defaultdict(set)
    enst_to_seq:  Dict[str, str]      = {}
    for header, seq in _iter_fasta(path):
        m = _RE_ENST.search(header)
        if not m or not seq:
            continue
        enst = m.group(1)
        seq_to_ensts[seq].add(enst)
        enst_to_seq[enst] = seq
    return dict(seq_to_ensts), enst_to_seq

def parse_refseq_pep(path: str):
    """Returns seq_to_nps: {seq -> set(NP.v)}, np_to_seq: {NP.v -> seq}"""
    seq_to_nps: Dict[str, Set[str]] = defaultdict(set)
    np_to_seq:  Dict[str, str]      = {}
    for header, seq in _iter_fasta(path):
        m = _RE_NP.match(header)
        if not m or not seq:
            continue
        np_ = m.group(1)
        seq_to_nps[seq].add(np_)
        np_to_seq[np_] = seq
    return dict(seq_to_nps), np_to_seq

# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------
def parse_gtf(path: str) -> Dict[str, Dict]:
    """Returns {ENST.v -> {gene_id, gene_symbol}} from transcript lines only."""
    enst_to_gene: Dict[str, Dict] = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if "\ttranscript\t" not in line:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            m_tx  = re.search(r'transcript_id "([^"]+)"', attrs)
            m_tv  = re.search(r'transcript_version "([^"]+)"', attrs)
            m_gid = re.search(r'gene_id "([^"]+)"', attrs)
            m_gn  = re.search(r'gene_name "([^"]+)"', attrs)
            if not m_tx:
                continue
            tid = m_tx.group(1)
            if m_tv:
                tid = f"{tid}.{m_tv.group(1)}"
            gid   = m_gid.group(1) if m_gid else ""
            gname = m_gn.group(1)  if m_gn  else ""
            gname = GENE_SYMBOL_UPDATES.get(gname, gname)
            enst_to_gene[tid] = {"gene_id": gid, "gene_symbol": gname}
    return enst_to_gene

# ---------------------------------------------------------------------------
# NP → NM mapping via gene2accession (25 GB — stream with cache)
# ---------------------------------------------------------------------------
def parse_refseq_tsv(path: str) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """
    Read pre-built refseq_transcripts_grch3[78].tsv (from build_gene_set.py).
    Columns: gene_symbol, entrez_id, hgnc_id, ensembl_gene_id,
             refseq_transcript_id, refseq_protein_id, protein_sequence
    Returns (nm_to_seq, gene_to_nms).
    """
    nm_to_seq: Dict[str, str] = {}
    gene_to_nms: Dict[str, List[str]] = defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        nm_col   = idx["refseq_transcript_id"]
        seq_col  = idx["protein_sequence"]
        gene_col = idx["gene_symbol"]
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(nm_col, seq_col, gene_col):
                continue
            nm   = fields[nm_col].strip()
            seq  = fields[seq_col].strip()
            gene = fields[gene_col].strip()
            gene = GENE_SYMBOL_UPDATES.get(gene, gene)
            if nm and seq:
                nm_to_seq[nm] = seq
            if nm and gene:
                gene_to_nms[gene].append(nm)
    # Deduplicate while preserving order
    for k in gene_to_nms:
        gene_to_nms[k] = sorted(set(gene_to_nms[k]), key=lambda x: (_nm_num(x), x))
    return nm_to_seq, dict(gene_to_nms)

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
        nm_col     = next((idx[k] for k in idx if "refseq" in k and "nuc" in k), None)
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
    Columns: Gene, MANE TYPE, Ensembl StableID GRCh38,
             RefSeq StableID GRCh38 / GRCh37, Ensembl StableID GRCh37 (Not MANE), ...
    Returns {gene -> {grch38_enst, grch37_enst, nm, mane_type}}
    """
    result: Dict[str, Dict] = {}
    with _open(path) as fh:
        header_line = fh.readline().rstrip("\n")
        sep = "\t" if "\t" in header_line else ","
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
                result[gene] = {"grch38_enst": grch38_enst, "grch37_enst": grch37_enst,
                                "nm": nm, "mane_type": mane_type}
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
    """Returns {gene -> isoform_override} — the first column of the Iv7 overrides file.
    Value can be a versioned ENST (e.g. ENST00000331340) or NM (e.g. NM_000038.5)."""
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
    Columns: gene_name  refseq_id  enst_id  note
    Returns {gene -> {refseq_id, enst_id, note}}
    """
    result: Dict[str, Dict] = {}
    with open(path, "r", encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            line = line.rstrip("\n")
            parts = line.split("\t")
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
    Columns: gene_name  refseq_id  enst_id  note
    Returns {gene -> {refseq_id, enst_id, note}}
    """
    result: Dict[str, Dict] = {}
    with open(path, "r", encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            line = line.rstrip("\n")
            parts = line.split("\t")
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

def parse_oncokb(path: str) -> Dict[str, Dict]:
    """
    Supports both old format (hugo_symbol, ...) and new format (with id, entrez_gene_id, etc.)
    Returns {gene -> {grch38_enst, grch38_nm}} for GRCh38 rows.
    """
    result: Dict[str, Dict] = {}
    with open(path, "r", encoding="utf-8") as fh:
        header_line = fh.readline().rstrip("\n")
        sep = "\t" if "\t" in header_line else ","
        cols = header_line.split(sep)
        for line in fh:
            fields = line.rstrip("\n").split(sep)
            row = dict(zip(cols, fields))
            # Support both column name variants; strip CSV quotes
            gene   = (row.get("hugo_symbol") or row.get("gene_symbol") or "").strip().strip('"')
            genome = (row.get("reference_genome") or "").strip().strip('"')
            enst   = (row.get("ensembl_transcript_id") or "").strip().strip('"')
            nm     = (row.get("reference_sequence_id") or "").strip().strip('"')
            # Normalize CDKN2A (p14) -> CDKN2A(p14) (remove space before paren)
            gene = re.sub(r'\s+\(', '(', gene)
            gene   = GENE_SYMBOL_UPDATES.get(gene, gene)
            if gene and genome == "GRCh38":
                result[gene] = {"grch38_enst": enst, "grch38_nm": nm}
    return result

def parse_oncokb_isoform(path: str) -> Dict[str, Dict]:
    """
    oncokb_isoform.tsv — GRCh37+38 transcripts + entrez gene ID.
    ENST IDs are unversioned.
    Returns {gene -> {entrez_gene_id, grch37_enst_base, grch37_nm, grch38_enst_base, grch38_nm, gene_type}}
    """
    result: Dict[str, Dict] = {}
    with open(path, "r", encoding="utf-8") as fh:
        header_line = fh.readline().rstrip("\n")
        cols = header_line.split("\t")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            row = dict(zip(cols, fields))
            gene = row.get("Hugo Symbol", "").strip()
            # Normalize CDKN2A (p14) -> CDKN2A(p14)
            gene = re.sub(r'\s+\(', '(', gene)
            gene = GENE_SYMBOL_UPDATES.get(gene, gene)
            entrez = row.get("Entrez Gene ID", "").strip()
            grch37_enst = row.get("GRCh37 Isoform", "").strip()
            grch37_nm = row.get("GRCh37 RefSeq", "").strip()
            grch38_enst = row.get("GRCh38 Isoform", "").strip()
            grch38_nm = row.get("GRCh38 RefSeq", "").strip()
            gene_type = row.get("Gene Type", "").strip()
            if gene:
                result[gene] = {
                    "entrez_gene_id": entrez,
                    "grch37_enst_base": grch37_enst,
                    "grch37_nm": grch37_nm,
                    "grch38_enst_base": grch38_enst,
                    "grch38_nm": grch38_nm,
                    "gene_type": gene_type,
                }
    return result

def parse_hgnc(path: str) -> Tuple[Dict[str, Dict], Dict[str, str], Dict[str, str]]:
    """
    Parse HGNC complete set TSV.
    Returns:
      hgnc_map:   {symbol -> {hgnc_id, entrez_id}}  (approved + prev)
      alias_map:  {alias_symbol -> canonical_symbol}
      prev_map:   {prev_symbol  -> canonical_symbol}
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
            # alias_symbol column
            if alias_col >= 0 and alias_col < len(fields):
                raw_alias = fields[alias_col].strip()
                if raw_alias and raw_alias != "nan":
                    for a in raw_alias.split("|"):
                        a = a.strip().strip('"')
                        if a and a not in hgnc_map and a not in alias_map:
                            alias_map[a] = symbol
            # prev_symbol column
            if prev_col < len(fields):
                raw_prev = fields[prev_col].strip()
                if raw_prev and raw_prev != "nan":
                    for ps in raw_prev.split("|"):
                        ps = ps.strip().strip('"')
                        if ps and ps not in hgnc_map and ps not in prev_map:
                            prev_map[ps] = symbol
                            # Also add to hgnc_map so existing callers still work
                            hgnc_map.setdefault(ps, info)
    return hgnc_map, alias_map, prev_map


def resolve_gene_symbol(gene: str, hgnc_map: Dict[str, Dict],
                        alias_map: Dict[str, str], prev_map: Dict[str, str]) -> str:
    """Resolve an outdated/alias gene name to its current approved HGNC symbol."""
    if gene in hgnc_map:
        return gene
    if gene in alias_map:
        return alias_map[gene]
    if gene in prev_map:
        return prev_map[gene]
    return GENE_SYMBOL_UPDATES.get(gene, gene)


def resolve_enst_version(enst_base: str, enst_to_seq: Dict[str, str]) -> str:
    """Resolve unversioned ENST to versioned ENST from FASTA/GTF data."""
    if not enst_base:
        return ""
    # If it already has a version and exists, return as-is
    if "." in enst_base and enst_base in enst_to_seq:
        return enst_base
    base = enst_base.split(".")[0]
    matches = [k for k in enst_to_seq if k.startswith(base + ".")]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        # Pick highest version
        return max(matches, key=_enst_ver)
    return enst_base

def parse_release221(path: str) -> Dict[str, Set[str]]:
    """
    release221.accession2geneid: taxid, entrez_gene_id, accession, protein_accession
    Stream-filter to taxid=9606, build {entrez_gene_id -> set(NM.v)}.
    Cached as pickle.
    """
    cache_path = os.path.join(CACHE_DIR, "release221_entrez_to_nms.pkl")
    if os.path.exists(cache_path):
        print("  Loading release221 entrez->NM map from cache…")
        with open(cache_path, "rb") as f:
            return pickle.load(f)
    print("  Streaming release221.accession2geneid (takes a few minutes)…")
    entrez_to_nms: Dict[str, Set[str]] = defaultdict(set)
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            if fields[0] != "9606":
                continue
            entrez = fields[1]
            accession = fields[2].strip()
            if accession.startswith("NM_"):
                entrez_to_nms[entrez].add(accession)
    result = dict(entrez_to_nms)
    with open(cache_path, "wb") as f:
        pickle.dump(result, f)
    print(f"  Cached {len(result):,} entrez->NM mappings")
    return result

# ---------------------------------------------------------------------------
# Similarity — uses difflib.SequenceMatcher for alignment-aware identity
# Correctly handles substitutions, insertions, and deletions.
# ---------------------------------------------------------------------------
def compute_similarity(seq1: str, seq2: str) -> Dict:
    if not seq1 or not seq2:
        return {"pct": None, "diff_count": 0}
    sm = difflib.SequenceMatcher(None, seq1, seq2, autojunk=False)
    pct = round(sm.ratio() * 100, 4)
    # Count individual non-equal operations
    diff_count = sum(
        max(i2 - i1, j2 - j1)
        for tag, i1, i2, j1, j2 in sm.get_opcodes()
        if tag != "equal"
    )
    return {"pct": pct, "diff_count": diff_count}

# ---------------------------------------------------------------------------
# ID helpers — "smallest" = lowest numeric value in the accession
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

# ---------------------------------------------------------------------------
# Best-match selection
# ---------------------------------------------------------------------------
def best_enst_for_gene(
    gene: str,
    ref_seq: str,
    gene_to_ensts: Dict[str, List[str]],
    enst_to_seq: Dict[str, str],
    mane_grch38: Dict[str, Dict],
    fallback_enst: str = "",
) -> Tuple[str, float, str]:
    """Returns (best_enst, pct_identity, source_label)."""
    if not ref_seq:
        return fallback_enst, 0.0, "No GRCh37 reference sequence available"

    candidates = gene_to_ensts.get(gene, [])
    if not candidates:
        return fallback_enst, 0.0, "No Ensembl GRCh38 transcripts found for gene"

    exact = [e for e in candidates if enst_to_seq.get(e, "") == ref_seq]
    if exact:
        # Priority: MANE Select > MANE Plus Clinical > same base > smallest ID
        for e in exact:
            if mane_grch38.get(e, {}).get("status") == "MANE Select":
                return e, 100.0, "MANE Select"
        for e in exact:
            if "MANE Plus Clinical" in mane_grch38.get(e, {}).get("status", ""):
                return e, 100.0, "MANE Plus Clinical"
        if fallback_enst:
            base = fallback_enst.split(".")[0]
            for e in exact:
                if e.split(".")[0] == base:
                    return e, 100.0, "Same ENST base ID (different version)"
        # Smallest numeric ENST ID
        best = min(exact, key=lambda e: (_enst_num(e), _enst_ver(e)))
        return best, 100.0, "Exact sequence match (smallest ID)"

    # No exact match — find most similar
    best_e, best_pct = fallback_enst, 0.0
    for e in candidates:
        seq = enst_to_seq.get(e, "")
        if not seq:
            continue
        sim = compute_similarity(ref_seq, seq)
        if sim["pct"] is not None and sim["pct"] > best_pct:
            best_pct = sim["pct"]
            best_e   = e
    label = f"Best available match ({best_pct:.1f}% identity)"
    return best_e, best_pct, label

def best_nm_for_gene(
    gene: str,
    ref_seq: str,
    gene_to_nms: Dict[str, List[str]],
    nm_to_seq: Dict[str, str],
    clinical_nm_base: Optional[str],
    is_germline: bool,
    mane_nm: str = "",
    oncokb_nm: str = "",
) -> Tuple[str, float, str]:
    """Returns (best_nm, pct_identity, source_label)."""
    if not ref_seq:
        return "", 0.0, "no reference sequence"

    candidates = gene_to_nms.get(gene, [])

    # Germline / clinical: look for the locked NM specifically
    if clinical_nm_base:
        locked = [n for n in candidates if n.split(".")[0] == clinical_nm_base]
        locked_exact = [n for n in locked if nm_to_seq.get(n, "") == ref_seq]
        if locked_exact:
            best = min(locked_exact, key=_nm_num)
            label = "Germline locked (Iv7)" if is_germline else "Clinical override (Iv7)"
            return best, 100.0, label
        if locked:
            best = min(locked, key=_nm_num)
            sim  = compute_similarity(ref_seq, nm_to_seq.get(best, ""))
            pct  = sim["pct"] if sim["pct"] is not None else 0.0
            label = f"Clinical NM (partial match {pct:.1f}%)"
            return best, pct, label

    # MANE NM exact match
    if mane_nm and candidates:
        mane_base = mane_nm.split(".")[0]
        mane_matches = [n for n in candidates if n.split(".")[0] == mane_base and nm_to_seq.get(n, "") == ref_seq]
        if mane_matches:
            return min(mane_matches, key=_nm_num), 100.0, "MANE Select"

    # OncoKB NM exact match
    if oncokb_nm and candidates:
        okb_base = oncokb_nm.split(".")[0]
        okb_matches = [n for n in candidates if n.split(".")[0] == okb_base and nm_to_seq.get(n, "") == ref_seq]
        if okb_matches:
            return min(okb_matches, key=_nm_num), 100.0, "OncoKB reference (exact match)"

    # Any exact match → smallest NM
    exact = [n for n in candidates if nm_to_seq.get(n, "") == ref_seq]
    if exact:
        return min(exact, key=_nm_num), 100.0, "Exact sequence match (smallest ID)"

    # Partial match → best then smallest
    best_n, best_pct = "", 0.0
    for n in candidates:
        seq = nm_to_seq.get(n, "")
        if not seq:
            continue
        sim = compute_similarity(ref_seq, seq)
        if sim["pct"] is not None and sim["pct"] > best_pct:
            best_pct = sim["pct"]
            best_n   = n
    label = f"Best available match ({best_pct:.1f}% identity)" if best_n else "not found"
    return best_n, best_pct, label

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=== Transcript Diff Pipeline ===")

    print("\n[1/11] Parsing Ensembl GRCh37 FASTA…")
    grch37_pep = finput("Homo_sapiens.GRCh37.pep.all.fa")
    if not os.path.exists(grch37_pep): grch37_pep += ".gz"
    seq_to_enst37, enst37_to_seq = parse_ensembl_pep(grch37_pep)
    print(f"  {len(enst37_to_seq):,} GRCh37 ENSTs")

    print("[2/11] Parsing Ensembl GRCh38 FASTA…")
    grch38_pep = finput("Homo_sapiens.GRCh38.pep.all.fa")
    if not os.path.exists(grch38_pep): grch38_pep += ".gz"
    seq_to_enst38, enst38_to_seq = parse_ensembl_pep(grch38_pep)
    print(f"  {len(enst38_to_seq):,} GRCh38 ENSTs")

    print("[3/11] Parsing RefSeq TSVs (from build_gene_set.py output)…")
    nm37_to_seq, gene_to_nm37_tsv = parse_refseq_tsv(
        os.path.join(ROOT, "files", "output", "refseq_transcripts_grch37.tsv")
    )
    nm38_to_seq, gene_to_nm38_tsv = parse_refseq_tsv(
        os.path.join(ROOT, "files", "output", "refseq_transcripts_grch38.tsv")
    )
    print(f"  {len(nm37_to_seq):,} GRCh37 NMs, {len(nm38_to_seq):,} GRCh38 NMs")

    print("[4/11] Parsing GTF files…")
    gtf37 = finput("Homo_sapiens.GRCh37.87.gtf")
    if not os.path.exists(gtf37): gtf37 += ".gz"
    enst37_to_gene = parse_gtf(gtf37)
    gtf38 = finput("Homo_sapiens.GRCh38.111.gtf")
    if not os.path.exists(gtf38): gtf38 += ".gz"
    enst38_to_gene = parse_gtf(gtf38)
    print(f"  {len(enst37_to_gene):,} GRCh37, {len(enst38_to_gene):,} GRCh38 transcripts")

    print("[5/11] Parsing MANE files…")
    mane38_path = fisoform("MANE.GRCh38.v1.2.summary.txt")
    if not os.path.exists(mane38_path): mane38_path += ".gz"
    mane_grch38 = parse_mane_grch38(mane38_path)
    mane_grch37_csv = parse_mane_grch37_csv(fisoform("MANE_GRCh37_list_filtered.csv"))
    print(f"  {len(mane_grch38):,} MANE GRCh38, {len(mane_grch37_csv):,} MANE GRCh37 entries")

    print("[6/11] Parsing clinical / annotation files…")
    iv7          = parse_iv7_overrides(fisoform("Iv7_dmp_isoform_merged_overrides.txt"))
    iv7_override = parse_iv7_isoform_override(fisoform("Iv7_dmp_isoform_merged_overrides.txt"))
    mskcc        = parse_mskcc_overrides(fisoform("isoform_overrides_at_mskcc_grch37.txt"))
    mskcc38      = parse_mskcc_grch38_overrides(fisoform("isoform_overrides_at_mskcc_grch38.txt"))
    germline     = parse_germline(fisoform("germline_panel_94.txt"))
    print(f"  Iv7: {len(iv7)}, MSKCC37: {len(mskcc)}, MSKCC38: {len(mskcc38)}, germline: {len(germline)}")

    print("[7/11] Parsing oncokb_isoform_versioned.tsv…")
    oncokb_isoform = parse_oncokb_isoform(fisoform("oncokb_isoform_versioned.tsv"))
    print(f"  OncoKB isoform entries: {len(oncokb_isoform)}")

    print("[8/11] Parsing HGNC complete set…")
    hgnc_path = fisoform("hgnc_complete_set_oct_07_2025.txt")
    hgnc_map, alias_map, prev_map = parse_hgnc(hgnc_path)
    print(f"  HGNC entries: {len(hgnc_map)}, aliases: {len(alias_map)}, prev: {len(prev_map)}")

    print("[9/11] Building gene -> transcript indexes…")
    # gene -> [ENST.v] for each assembly
    gene_to_enst37: Dict[str, List[str]] = defaultdict(list)
    for enst, info in enst37_to_gene.items():
        g = info["gene_symbol"]
        if g:
            gene_to_enst37[g].append(enst)

    gene_to_enst38: Dict[str, List[str]] = defaultdict(list)
    for enst, info in enst38_to_gene.items():
        g = info["gene_symbol"]
        if g:
            gene_to_enst38[g].append(enst)

    # gene -> [NM.v] from pre-built RefSeq TSVs (GFF3-based, comprehensive)
    gene_to_nm37: Dict[str, List[str]] = defaultdict(list, gene_to_nm37_tsv)
    gene_to_nm38: Dict[str, List[str]] = defaultdict(list, gene_to_nm38_tsv)

    # Remove duplicates while preserving lists
    for d in (gene_to_nm37, gene_to_nm38, gene_to_enst37, gene_to_enst38):
        for k in d:
            d[k] = sorted(set(d[k]), key=lambda x: (_nm_num(x) if x.startswith("NM") else _enst_num(x), x))

    all_genes: Set[str] = set(iv7) | set(mskcc) | set(mskcc38) | set(oncokb_isoform)
    print(f"  Gene scope: {len(all_genes)}")

    missing_in_mskcc = set(iv7) - set(mskcc)
    if missing_in_mskcc:
        print(f"  WARNING: {len(missing_in_mskcc)} Iv7 genes not in MSKCC: {sorted(missing_in_mskcc)}")

    print("[10/11] Building gene records…")

    genes_out = []
    warnings  = []

    for gene in sorted(all_genes):
        mskcc_entry   = mskcc.get(gene, {})
        mskcc38_entry = mskcc38.get(gene, {})
        iv7_nm_base   = iv7.get(gene)
        is_germline   = gene in germline
        is_clinical   = gene in iv7
        is_oncokb     = gene in oncokb_isoform
        okb_iso       = oncokb_isoform.get(gene, {})

        # ---- HGNC + Entrez Gene ID (HGNC first, fallback oncokb_isoform) ----
        hgnc_entry = hgnc_map.get(gene, {})
        hgnc_id = hgnc_entry.get("hgnc_id", "")
        entrez_gene_id = hgnc_entry.get("entrez_id", "") or okb_iso.get("entrez_gene_id", "")

        # ---- GRCh37 ENST (from MSKCC curated) ----
        grch37_enst = mskcc_entry.get("enst_id", "")
        # Resolve to exact versioned ID if needed
        if grch37_enst and grch37_enst not in enst37_to_seq:
            base = grch37_enst.split(".")[0]
            matched = next((k for k in enst37_to_seq if k.startswith(base + ".")), None)
            grch37_enst = matched or grch37_enst
        grch37_enst_sources = []
        if mskcc_entry:
            grch37_enst_sources.append("Curated (MSKCC isoform overrides)")

        # Fallback: MANE GRCh37
        if not grch37_enst:
            mane37 = mane_grch37_csv.get(gene, {})
            grch37_enst = mane37.get("grch37_enst", "")
            if grch37_enst:
                grch37_enst_sources = ["MANE (GRCh37 lifted)"]
        else:
            # Check if also in MANE GRCh37
            mane37 = mane_grch37_csv.get(gene, {})
            if mane37.get("grch37_enst", "") and grch37_enst.split(".")[0] == mane37.get("grch37_enst", "").split(".")[0]:
                grch37_enst_sources.append("MANE GRCh37 (lifted)")

        # Check if also in oncokb_isoform GRCh37
        if grch37_enst and okb_iso.get("grch37_enst_base"):
            if grch37_enst.split(".")[0] == okb_iso["grch37_enst_base"].split(".")[0]:
                grch37_enst_sources.append("OncoKB isoform")

        grch37_enst_source = "; ".join(grch37_enst_sources) if grch37_enst_sources else ""
        grch37_seq = enst37_to_seq.get(grch37_enst, "")

        # ---- GRCh37 NM ----
        grch37_nm        = mskcc_entry.get("refseq_id", "")
        grch37_nm_sources = []
        if grch37_nm:
            grch37_nm_sources.append("Curated (MSKCC isoform overrides)")
            # Check if also matches Iv7 clinical
            if iv7_nm_base and grch37_nm.split(".")[0] == iv7_nm_base:
                grch37_nm_sources.append("Iv7 clinical override")
            # Try to resolve to a versioned NM that exists in our seq map
            if grch37_nm not in nm37_to_seq:
                base = grch37_nm.split(".")[0]
                candidates_nm37 = [n for n in nm37_to_seq if n.split(".")[0] == base]
                if candidates_nm37:
                    exact_matches = [n for n in candidates_nm37 if grch37_seq and nm37_to_seq[n] == grch37_seq]
                    if exact_matches:
                        grch37_nm = min(exact_matches, key=_nm_num)
                    else:
                        grch37_nm = min(candidates_nm37, key=_nm_num)
        else:
            # Fallback: sequence-based search
            grch37_nm, _, grch37_nm_label = best_nm_for_gene(
                gene, grch37_seq, gene_to_nm37, nm37_to_seq,
                iv7_nm_base, is_germline
            )
            if grch37_nm_label:
                grch37_nm_sources.append(grch37_nm_label)
        grch37_nm_source = "; ".join(grch37_nm_sources) if grch37_nm_sources else ""

        # ---- GRCh38 ENST ----
        grch38_enst_hint = okb_iso.get("grch38_enst_base", "")
        if grch38_enst_hint and grch38_enst_hint not in enst38_to_seq:
            base = grch38_enst_hint.split(".")[0]
            grch38_enst_hint = next((k for k in enst38_to_seq if k.startswith(base + ".")), "")

        grch38_enst, grch38_enst_pct, grch38_enst_label = best_enst_for_gene(
            gene, grch37_seq, gene_to_enst38, enst38_to_seq, mane_grch38, grch38_enst_hint
        )
        if not grch38_enst and grch38_enst_hint:
            grch38_enst, grch38_enst_label = grch38_enst_hint, "OncoKB reference (no seq match)"

        # Build multi-source label for GRCh38 ENST
        grch38_enst_sources = [grch38_enst_label] if grch38_enst_label else []
        if grch38_enst and okb_iso.get("grch38_enst_base"):
            if grch38_enst.split(".")[0] == okb_iso["grch38_enst_base"].split(".")[0]:
                if "OncoKB" not in "; ".join(grch38_enst_sources):
                    grch38_enst_sources.append("OncoKB isoform")
        grch38_enst_source = "; ".join(grch38_enst_sources)

        # ---- GRCh38 NM ----
        mane_entry  = mane_grch38.get(grch38_enst, {}) if grch38_enst else {}
        mane38_nm   = mane_entry.get("nm", "")
        oncokb_nm38 = okb_iso.get("grch38_nm", "")

        grch38_nm, grch38_nm_pct, grch38_nm_label = best_nm_for_gene(
            gene, grch37_seq, gene_to_nm38, nm38_to_seq,
            iv7_nm_base, is_germline, mane38_nm, oncokb_nm38
        )
        # Build multi-source label for GRCh38 NM
        grch38_nm_sources = [grch38_nm_label] if grch38_nm_label else []
        if grch38_nm and iv7_nm_base and grch38_nm.split(".")[0] == iv7_nm_base:
            if "clinical" not in grch38_nm_label.lower() and "Iv7" not in grch38_nm_label:
                grch38_nm_sources.append("Iv7 clinical override")
        grch38_nm_source = "; ".join(grch38_nm_sources)

        # ---- MANE GRCh38 ----
        mane38_status = mane_entry.get("status", "")
        mane38_enst   = grch38_enst if "MANE" in mane38_status else ""
        if not mane38_enst:
            for e in gene_to_enst38.get(gene, []):
                s = mane_grch38.get(e, {}).get("status", "")
                if "MANE Select" in s:
                    mane38_enst  = e
                    mane38_nm    = mane_grch38[e].get("nm", "")
                    mane38_status = s
                    break

        # ---- MANE GRCh37 validated ----
        mane37_entry       = mane_grch37_csv.get(gene, {})
        mane37_grch38_enst = mane37_entry.get("grch38_enst", "")
        mane37_grch37_enst = mane37_entry.get("grch37_enst", "")
        mane_grch37_valid  = mane37_grch38_enst in mane_grch38 if mane37_grch38_enst else False

        curated_note = mskcc_entry.get("note", "")

        # ---- MSKCC GRCh38 ENST / NM ----
        mskcc38_enst_raw = mskcc38_entry.get("enst_id", "")
        mskcc38_nm_raw   = mskcc38_entry.get("refseq_id", "")
        # Resolve versioned ENST
        if mskcc38_enst_raw and mskcc38_enst_raw not in enst38_to_seq:
            base38 = mskcc38_enst_raw.split(".")[0]
            matched38 = next((k for k in enst38_to_seq if k.startswith(base38 + ".")), None)
            mskcc38_enst_raw = matched38 or mskcc38_enst_raw
        # Resolve versioned NM
        if mskcc38_nm_raw and mskcc38_nm_raw not in nm38_to_seq:
            base38nm = mskcc38_nm_raw.split(".")[0]
            candidates38nm = [n for n in nm38_to_seq if n.split(".")[0] == base38nm]
            if candidates38nm:
                mskcc38_nm_raw = min(candidates38nm, key=_nm_num)

        # ---- mane_only flag ----
        # True when the gene has no curated entry except MANE
        _has_curated = bool(
            mskcc_entry or mskcc38_entry or iv7_nm_base or
            okb_iso.get("grch37_enst_base") or okb_iso.get("grch38_enst_base")
        )
        mane_only = not _has_curated

        # ---- OncoKB best_match (resolved versioned IDs from oncokb_isoform.tsv) ----
        oncokb_best = None
        if is_oncokb:
            okb37_enst = resolve_enst_version(okb_iso.get("grch37_enst_base", ""), enst37_to_seq)
            okb38_enst = resolve_enst_version(okb_iso.get("grch38_enst_base", ""), enst38_to_seq)
            okb37_nm = okb_iso.get("grch37_nm", "")
            okb38_nm = okb_iso.get("grch38_nm", "")
            # Resolve NM versions if needed
            def _resolve_nm(nm_val: str, nm_seq_map: Dict[str, str]) -> str:
                if nm_val and nm_val not in nm_seq_map:
                    base = nm_val.split(".")[0]
                    matches = [k for k in nm_seq_map if k.startswith(base + ".")]
                    if matches:
                        return max(matches, key=_enst_ver)
                return nm_val
            okb37_nm = _resolve_nm(okb37_nm, nm37_to_seq)
            okb38_nm = _resolve_nm(okb38_nm, nm38_to_seq)
            oncokb_best = {
                "grch37_enst": okb37_enst or "",
                "grch38_enst": okb38_enst or "",
                "grch37_nm": okb37_nm or "",
                "grch38_nm": okb38_nm or "",
            }

        # ---- Warnings ----
        if is_germline and grch37_seq:
            s38 = enst38_to_seq.get(grch38_enst, "")
            if s38 and s38 != grch37_seq:
                warnings.append(f"GERMLINE_MISMATCH {gene}: GRCh37 ENST seq != GRCh38 ENST seq")
        if is_clinical and not grch37_seq:
            warnings.append(f"NO_SEQ {gene}: clinical gene has no GRCh37 ENST sequence")
        if grch38_enst_pct < 100 and grch38_enst_pct > 0:
            warnings.append(f"PARTIAL_MATCH {gene}: GRCh38 ENST {grch38_enst_pct:.1f}% identity")
        if is_clinical and not grch37_nm:
            warnings.append(f"NO_GRCh37_NM {gene}: clinical gene has no GRCh37 NM")
        if is_clinical and not grch38_nm:
            warnings.append(f"NO_GRCh38_NM {gene}: clinical gene has no GRCh38 NM")

        # ---- Collect all transcripts with source labels ----
        all_transcripts = []
        seen_keys: Set[Tuple[str, str]] = set()  # (tid, assembly) to allow same ID across assemblies

        def add_transcript(tid: str, t_type: str, assembly: str,
                           seq_map: Dict[str, str], source: str, is_primary: bool = False):
            key = (tid, assembly)
            if not tid or key in seen_keys:
                return
            seen_keys.add(key)
            seq = seq_map.get(tid, "")
            if not seq:
                base = tid.split(".")[0]
                for k, v in seq_map.items():
                    if k.startswith(base + "."):
                        seq = v
                        break
            all_transcripts.append({
                "id":         tid,
                "type":       t_type,
                "assembly":   assembly,
                "sequence":   seq,
                "length":     len(seq),
                "source":     source,
                "is_primary": is_primary,
            })

        # Primary transcripts (is_primary = True)
        add_transcript(grch37_enst, "ensembl", "GRCh37", enst37_to_seq, grch37_enst_source or "Curated (MSKCC)", True)
        add_transcript(grch38_enst, "ensembl", "GRCh38", enst38_to_seq, grch38_enst_source, True)
        add_transcript(grch37_nm,   "refseq",  "GRCh37", nm37_to_seq,   grch37_nm_source,   True)
        add_transcript(grch38_nm,   "refseq",  "GRCh38", nm38_to_seq,   grch38_nm_source,   True)
        if mane38_enst and mane38_enst != grch38_enst:
            add_transcript(mane38_enst, "ensembl", "GRCh38", enst38_to_seq, "MANE Select (reference)", False)
        # MSKCC GRCh38 isoform transcripts
        add_transcript(mskcc38_enst_raw, "ensembl", "GRCh38", enst38_to_seq, "MSKCC isoform override (GRCh38)")
        add_transcript(mskcc38_nm_raw,   "refseq",  "GRCh38", nm38_to_seq,   "MSKCC isoform override (GRCh38)")

        # OncoKB isoform transcripts (if different from primaries)
        if oncokb_best:
            add_transcript(oncokb_best["grch37_enst"], "ensembl", "GRCh37", enst37_to_seq, "OncoKB isoform (GRCh37)")
            add_transcript(oncokb_best["grch38_enst"], "ensembl", "GRCh38", enst38_to_seq, "OncoKB isoform (GRCh38)")
            add_transcript(oncokb_best["grch37_nm"],   "refseq",  "GRCh37", nm37_to_seq,   "OncoKB isoform (GRCh37)")
            add_transcript(oncokb_best["grch38_nm"],   "refseq",  "GRCh38", nm38_to_seq,   "OncoKB isoform (GRCh38)")

        # All alternative Ensembl transcripts from GTF
        for e in gene_to_enst37.get(gene, []):
            add_transcript(e, "ensembl", "GRCh37", enst37_to_seq, "Alternative Ensembl GRCh37 transcript")
        for e in gene_to_enst38.get(gene, []):
            add_transcript(e, "ensembl", "GRCh38", enst38_to_seq, "Alternative Ensembl GRCh38 transcript")
        # All NMs from RefSeq TSVs (GFF3-based, comprehensive)
        for n in gene_to_nm37.get(gene, []):
            add_transcript(n, "refseq", "GRCh37", nm37_to_seq, "RefSeq GRCh37 NM (GFF3)")
        for n in gene_to_nm38.get(gene, []):
            add_transcript(n, "refseq", "GRCh38", nm38_to_seq, "RefSeq GRCh38 NM (GFF3)")

        # ---- Similarities (alignment-aware via SequenceMatcher) ----
        s37e = enst37_to_seq.get(grch37_enst, "")
        s38e = enst38_to_seq.get(grch38_enst, "") if grch38_enst else ""
        s37n = nm37_to_seq.get(grch37_nm, "")     if grch37_nm   else ""
        s38n = nm38_to_seq.get(grch38_nm, "")     if grch38_nm   else ""
        s_mane = enst38_to_seq.get(mane38_enst, "") if mane38_enst else ""

        similarities = {
            "grch37_enst_vs_grch37_refseq": compute_similarity(s37e, s37n),
            "grch38_enst_vs_grch38_refseq": compute_similarity(s38e, s38n),
            "grch37_enst_vs_grch38_enst":   compute_similarity(s37e, s38e),
            "grch37_enst_vs_mane":          compute_similarity(s37e, s_mane),
            "grch38_enst_vs_mane":          compute_similarity(s38e, s_mane),
        }

        gene_id = ""
        if grch37_enst and grch37_enst in enst37_to_gene:
            gene_id = enst37_to_gene[grch37_enst].get("gene_id", "")
        if not gene_id and grch38_enst and grch38_enst in enst38_to_gene:
            gene_id = enst38_to_gene[grch38_enst].get("gene_id", "")

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
            "mskcc38_nm":      mskcc38_nm_raw   or "",
            "mane_grch38": {
                "enst":   mane38_enst   or "",
                "nm":     mane38_nm     or "",
                "status": mane38_status or "",
            },
            "mane_grch37_enst": mane37_grch37_enst if mane_grch37_valid else "",
            "best_match": {
                "grch37_enst": grch37_enst or "",
                "grch38_enst": grch38_enst or "",
                "grch37_nm":   grch37_nm   or "",
                "grch38_nm":   grch38_nm   or "",
            },
            "transcripts":  all_transcripts,
            "similarities": similarities,
        }
        if oncokb_best:
            gene_record["oncokb_best_match"] = oncokb_best
        genes_out.append(gene_record)

    # ---- Collections ----
    print("[11/11] Building collections…")
    collections: Dict = {
        "oncokb": {
            g: {"grch38_enst": v.get("grch38_enst_base", ""), "grch38_nm": v.get("grch38_nm", "")}
            for g, v in oncokb_isoform.items()
        },
        "oncokb_isoform": {
            g: {
                "entrez_gene_id": v.get("entrez_gene_id", ""),
                "grch37_enst": v.get("grch37_enst_base", ""),
                "grch37_nm": v.get("grch37_nm", ""),
                "grch38_enst": v.get("grch38_enst_base", ""),
                "grch38_nm": v.get("grch38_nm", ""),
                "gene_type": v.get("gene_type", ""),
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
    for enst, info in mane_grch38.items():
        gene_info = enst38_to_gene.get(enst, {})
        g = gene_info.get("gene_symbol", "")
        if g:
            collections["mane"][g] = {
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
            print(f"  … and {len(warnings)-60} more")

if __name__ == "__main__":
    main()
