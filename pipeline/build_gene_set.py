#!/usr/bin/env python3
"""
Build gene set + transcript tables for Ensembl and RefSeq (GRCh37 and GRCh38).

Run from the project root:
    python3 pipeline/build_gene_set.py

Outputs (in files/output/):
    gene_set.tsv
    ensembl_transcripts_grch37.tsv
    ensembl_transcripts_grch38.tsv
    refseq_transcripts_grch37.tsv
    refseq_transcripts_grch38.tsv
"""

import csv
import os
import re
import sys
from collections import defaultdict
from typing import Dict, Optional, Set, Tuple

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES  = os.path.join(ROOT, "files")
OUTPUT = os.path.join(FILES, "output")
os.makedirs(OUTPUT, exist_ok=True)

def fp(*parts):
    return os.path.join(FILES, *parts)

def finput(*parts):
    return os.path.join(FILES, "input", *parts)

def fisoform(*parts):
    return os.path.join(FILES, "isoform", *parts)

def outp(*parts):
    return os.path.join(OUTPUT, *parts)

# ---------------------------------------------------------------------------
# Gene symbol normalisation (same as build_data.py)
# ---------------------------------------------------------------------------
GENE_SYMBOL_UPDATES = {
    "EIF2C1": "AGO1",
    "FTSJD1": "CMTR2",
    "NUT":    "NUTM1",
    "PAK7":   "PAK5",
    "WHSC1":  "NSD2",
    "STK19":  "WHR1",
}

def normalize_symbol(sym: str) -> str:
    sym = sym.strip()
    return GENE_SYMBOL_UPDATES.get(sym, sym)

# ---------------------------------------------------------------------------
# HGNC lookup table
# ---------------------------------------------------------------------------
def load_hgnc(path: str):
    """
    Returns:
        symbol_to_info:  {symbol -> {hgnc_id, entrez_id, ensembl_gene_id}}
        alias_to_symbols: {alias -> [approved_symbol, ...]}
        prev_to_symbols:  {prev_sym -> [approved_symbol, ...]}
        hgnc_id_to_info: {'HGNC:1100' -> {symbol, entrez_id, ensembl_gene_id}}
    """
    symbol_to_info = {}
    alias_to_symbols = defaultdict(list)
    prev_to_symbols  = defaultdict(list)
    hgnc_id_to_info  = {}  # 'HGNC:NNN' -> {symbol, entrez_id, ensembl_gene_id}

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym     = row.get("symbol", "").strip()
            hgnc_id = row.get("hgnc_id", "").strip()
            entrez  = row.get("entrez_id", "").strip()
            ensg    = row.get("ensembl_gene_id", "").strip()
            if not sym:
                continue
            # HGNC TSV sometimes stores entrez_id as float (e.g. '672.0') — normalise to int string
            if entrez and entrez.endswith(".0"):
                entrez = entrez[:-2]
            info = {
                "hgnc_id":         hgnc_id,
                "entrez_id":       entrez,
                "ensembl_gene_id": ensg,
            }
            symbol_to_info[sym] = info
            if hgnc_id:
                hgnc_id_to_info[hgnc_id] = {"symbol": sym, **info}

            # alias_symbol column: pipe-separated (or nan)
            for alias in _split_pipe(row.get("alias_symbol", "")):
                alias_to_symbols[alias].append(sym)
            # prev_symbol column: pipe-separated
            for prev in _split_pipe(row.get("prev_symbol", "")):
                prev_to_symbols[prev].append(sym)

    return symbol_to_info, alias_to_symbols, prev_to_symbols, hgnc_id_to_info


def _split_pipe(val: str):
    val = val.strip()
    if not val or val.lower() == "nan":
        return []
    return [x.strip() for x in val.split("|") if x.strip()]


def hgnc_lookup(symbol: str, symbol_to_info, alias_to_symbols, prev_to_symbols):
    """Return info dict or None."""
    if symbol in symbol_to_info:
        return symbol_to_info[symbol]
    # Try alias
    for approved in alias_to_symbols.get(symbol, []):
        if approved in symbol_to_info:
            return symbol_to_info[approved]
    # Try prev_symbol
    for approved in prev_to_symbols.get(symbol, []):
        if approved in symbol_to_info:
            return symbol_to_info[approved]
    return None

# ---------------------------------------------------------------------------
# Step 1 — Collect gene symbols from all 4 sources
# ---------------------------------------------------------------------------
def collect_gene_symbols() -> Set[str]:
    symbols: Set[str] = set()

    # 1a. Iv7_dmp_isoform_merged_overrides.txt — columns: #isoform_override, gene_name, dmp_refseq_id
    iv7_path = fisoform("Iv7_dmp_isoform_merged_overrides.txt")
    with open(iv7_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("gene_name", "").strip()
            if sym:
                symbols.add(normalize_symbol(sym))

    # 1b. oncokb_isoform_versioned.tsv — column: Hugo Symbol
    oncokb_path = fisoform("oncokb_isoform_versioned.tsv")
    with open(oncokb_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("Hugo Symbol", "").strip().strip('"')
            if sym:
                symbols.add(normalize_symbol(sym))

    # 1c. isoform_overrides_at_mskcc_grch37.txt — column: gene_name
    mskcc_path = fisoform("isoform_overrides_at_mskcc_grch37.txt")
    with open(mskcc_path, newline="", encoding="utf-8") as f:
        # File may have Windows line endings
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("gene_name", "").strip()
            if sym:
                symbols.add(normalize_symbol(sym))

    # 1d. MANE.GRCh38.v1.2.summary.txt — column: symbol
    mane_path = fisoform("MANE.GRCh38.v1.2.summary.txt")
    with open(mane_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("symbol", "").strip()
            if sym:
                symbols.add(normalize_symbol(sym))

    return symbols


# ---------------------------------------------------------------------------
# Step 2 — Build gene set TSV
# ---------------------------------------------------------------------------
def load_mane_fallback(mane_path: str, hgnc_id_to_info: Dict[str, dict] = None) -> Dict[str, dict]:
    """
    Build a {mane_symbol -> {entrez_id, ensembl_gene_id, hgnc_id, official_symbol}} index.
    MANE columns: #NCBI_GeneID  Ensembl_Gene  HGNC_ID  symbol ...
    If hgnc_id_to_info is provided, resolves the latest official gene symbol via HGNC_ID.
    """
    if hgnc_id_to_info is None:
        hgnc_id_to_info = {}

    fallback: Dict[str, dict] = {}
    with open(mane_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = row.get("symbol", "").strip()
            if not sym:
                continue
            # NCBI_GeneID format: 'GeneID:672'
            ncbi_raw = row.get("#NCBI_GeneID", "").strip()
            entrez = ncbi_raw.split(":", 1)[-1] if ":" in ncbi_raw else ncbi_raw
            ensg_raw = row.get("Ensembl_Gene", "").strip()
            # Ensembl_Gene may include version like 'ENSG00000012048.23' — strip version
            ensg = ensg_raw.split(".")[0] if ensg_raw else ""
            hgnc_raw = row.get("HGNC_ID", "").strip()  # already like 'HGNC:1100'

            # Resolve official symbol from HGNC if possible
            hgnc_entry = hgnc_id_to_info.get(hgnc_raw, {})
            official_symbol = hgnc_entry.get("symbol", "")  # latest approved symbol

            if sym not in fallback:  # keep first occurrence per mane symbol
                fallback[sym] = {
                    "entrez_id":       entrez,
                    "ensembl_gene_id": ensg,
                    "hgnc_id":         hgnc_raw,
                    "official_symbol": official_symbol,  # may be empty if HGNC not found
                }
    return fallback


def build_gene_set(symbols: Set[str], symbol_to_info, alias_to_symbols, prev_to_symbols,
                   mane_fallback: Dict[str, dict] = None):
    """
    Returns list of dicts with gene_symbol, entrez_id, hgnc_id, ensembl_gene_id.

    Fallback order for unresolved symbols:
      1. CDKN2A(p14) → use CDKN2A IDs, keep original symbol
      2. Any symbol found in MANE → use MANE NCBI_GeneID / Ensembl_Gene / HGNC_ID
      3. Truly unresolved → skip with warning
    """
    if mane_fallback is None:
        mane_fallback = {}

    gene_set = []
    not_found = []

    for sym in sorted(symbols):
        info = hgnc_lookup(sym, symbol_to_info, alias_to_symbols, prev_to_symbols)
        if info is not None:
            gene_set.append({
                "gene_symbol":     sym,
                "entrez_id":       info["entrez_id"],
                "hgnc_id":         info["hgnc_id"],
                "ensembl_gene_id": info["ensembl_gene_id"],
            })
            continue

        # Fallback 1: CDKN2A(p14) — use CDKN2A IDs
        if sym == "CDKN2A(p14)":
            cdkn2a_info = hgnc_lookup("CDKN2A", symbol_to_info, alias_to_symbols, prev_to_symbols)
            if cdkn2a_info:
                gene_set.append({
                    "gene_symbol":     sym,   # keep CDKN2A(p14)
                    "entrez_id":       cdkn2a_info["entrez_id"],
                    "hgnc_id":         cdkn2a_info["hgnc_id"],
                    "ensembl_gene_id": cdkn2a_info["ensembl_gene_id"],
                })
                print(f"  ℹ️  {sym}: used CDKN2A gene IDs")
                continue

        # Fallback 2: look up in MANE by symbol
        mane_info = mane_fallback.get(sym)
        if mane_info:
            official_sym = mane_info.get("official_symbol") or sym
            if official_sym != sym:
                print(f"  ℹ️  {sym}: resolved via MANE → official symbol '{official_sym}' "
                      f"(entrez={mane_info['entrez_id']})")
            else:
                print(f"  ℹ️  {sym}: resolved via MANE (entrez={mane_info['entrez_id']})")
            gene_set.append({
                "gene_symbol":     official_sym,
                "entrez_id":       mane_info["entrez_id"],
                "hgnc_id":         mane_info["hgnc_id"],
                "ensembl_gene_id": mane_info["ensembl_gene_id"],
            })
            continue

        not_found.append(sym)

    if not_found:
        print(f"\n⚠️  WARNING — {len(not_found)} gene(s) NOT found in HGNC or MANE and were SKIPPED:")
        for g in sorted(not_found):
            print(f"   {g}")
        print()

    return gene_set


def write_gene_set(gene_set, path):
    cols = ["gene_symbol", "entrez_id", "hgnc_id", "ensembl_gene_id"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(gene_set)
    print(f"  Wrote {len(gene_set)} genes → {path}")


# ---------------------------------------------------------------------------
# Step 3 — Ensembl transcript tables
# ---------------------------------------------------------------------------

def parse_ensembl_pep(fasta_path: str) -> Dict[str, str]:
    """
    Returns {ensp_versioned -> protein_sequence}
    Header example:
      >ENSP00000263100.2 pep:known ... transcript:ENST00000263100.8 gene:ENSG00000121410.12 ...
    """
    ensp_to_seq = {}
    current_id = None
    current_seq_parts = []

    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    ensp_to_seq[current_id] = "".join(current_seq_parts)
                # Parse ID — first token after ">"
                current_id = line[1:].split()[0]  # e.g. ENSP00000263100.2
                current_seq_parts = []
            else:
                if current_id is not None:
                    current_seq_parts.append(line)
        if current_id is not None:
            ensp_to_seq[current_id] = "".join(current_seq_parts)

    return ensp_to_seq


def parse_ensembl_gtf(gtf_path: str):
    """
    Parses CDS lines from an Ensembl GTF.
    Returns:
        ensg_to_transcripts: {ensg_id -> set of (enst_versioned, ensp_versioned)}
        ensg_to_symbol: {ensg_id -> gene_name}
        symbol_to_ensg: {gene_name -> ensg_id}  (first ENSG seen per symbol)
    """
    RE_ATTR = re.compile(r'(\w+) "([^"]+)"')

    ensg_to_transcripts: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)
    ensg_to_symbol: Dict[str, str] = {}
    symbol_to_ensg: Dict[str, str] = {}

    with open(gtf_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature != "CDS":
                continue

            attrs_str = parts[8]
            attrs = dict(RE_ATTR.findall(attrs_str))

            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", "")
            transcript_id = attrs.get("transcript_id", "")
            transcript_version = attrs.get("transcript_version", "")
            protein_id = attrs.get("protein_id", "")
            protein_version = attrs.get("protein_version", "")

            if not (gene_id and transcript_id and protein_id):
                continue

            enst_v = f"{transcript_id}.{transcript_version}" if transcript_version else transcript_id
            ensp_v = f"{protein_id}.{protein_version}" if protein_version else protein_id

            ensg_to_transcripts[gene_id].add((enst_v, ensp_v))
            if gene_name:
                ensg_to_symbol[gene_id] = gene_name
                if gene_name not in symbol_to_ensg:
                    symbol_to_ensg[gene_name] = gene_id

    return ensg_to_transcripts, ensg_to_symbol, symbol_to_ensg


def build_ensembl_table(gene_set, gtf_path: str, fasta_path: str, assembly: str,
                        sym_to_prev: Dict[str, list] = None):
    """
    Returns list of row dicts for the Ensembl transcript table.

    sym_to_prev: {current_symbol -> [prev_symbol, ...]} from HGNC.  When a
    gene's HGNC ENSG is absent from the GTF (common for genes renamed between
    GRCh37 and GRCh38), the function falls back to looking up the gene by its
    current symbol and then by each previous symbol in the GTF's gene-name
    index.  This recovers transcripts for genes like CCNQ (was FAM58A),
    H3C1 (was HIST1H3A), IKBKE (same name, different ENSG), etc.
    """
    print(f"  Parsing Ensembl GTF: {os.path.basename(gtf_path)} ...")
    ensg_to_transcripts, ensg_to_symbol, symbol_to_ensg = parse_ensembl_gtf(gtf_path)

    print(f"  Parsing Ensembl PEP FASTA: {os.path.basename(fasta_path)} ...")
    ensp_to_seq = parse_ensembl_pep(fasta_path)

    rows = []
    no_seq_count = 0
    fallback_count = 0

    for gene in gene_set:
        ensg = gene["ensembl_gene_id"]
        if not ensg:
            continue
        transcripts = ensg_to_transcripts.get(ensg, set())

        # Fallback: HGNC ENSG not present in this GTF (gene renamed or ENSG
        # reassigned between GRCh37 and GRCh38).  Try the GTF's own symbol
        # index using the current symbol and then each previous HGNC symbol.
        if not transcripts and sym_to_prev is not None:
            sym = gene["gene_symbol"]
            candidates = [sym] + list(sym_to_prev.get(sym, []))
            for candidate in candidates:
                fallback_ensg = symbol_to_ensg.get(candidate)
                if fallback_ensg and fallback_ensg != ensg:
                    transcripts = ensg_to_transcripts.get(fallback_ensg, set())
                    if transcripts:
                        gtf_name = ensg_to_symbol.get(fallback_ensg, candidate)
                        print(f"  ℹ️  {sym}: HGNC ENSG {ensg} not in GTF; "
                              f"resolved via GTF symbol '{gtf_name}' → {fallback_ensg} "
                              f"({len(transcripts)} transcript(s))")
                        fallback_count += 1
                        break

        for (enst_v, ensp_v) in sorted(transcripts):
            seq = ensp_to_seq.get(ensp_v, "")
            if not seq:
                no_seq_count += 1
            rows.append({
                "gene_symbol":          gene["gene_symbol"],
                "entrez_id":            gene["entrez_id"],
                "hgnc_id":              gene["hgnc_id"],
                "ensembl_gene_id":      ensg,
                "ensembl_transcript_id": enst_v,
                "ensembl_protein_id":   ensp_v,
                "protein_sequence":     seq,
            })

    if fallback_count:
        print(f"  ℹ️  {fallback_count} gene(s) resolved via GTF symbol fallback")
    if no_seq_count:
        print(f"  ⚠️  {no_seq_count} transcripts had no protein sequence in FASTA (may be non-coding)")

    return rows


def write_ensembl_table(rows, path):
    cols = [
        "gene_symbol", "entrez_id", "hgnc_id", "ensembl_gene_id",
        "ensembl_transcript_id", "ensembl_protein_id", "protein_sequence"
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows)} rows → {path}")


# ---------------------------------------------------------------------------
# Step 4 — RefSeq transcript tables
# ---------------------------------------------------------------------------

def _parse_gff_attrs(attrs_str: str) -> Dict[str, str]:
    """Parse GFF3 attribute string 'key=val;key=val' into dict."""
    result = {}
    for part in attrs_str.split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            result[k.strip()] = v.strip()
    return result


def _parse_dbxref(dbxref_str: str) -> Dict[str, str]:
    """Parse Dbxref like 'GeneID:1234,GenBank:NM_xxx.1,HGNC:HGNC:5' into dict."""
    result = {}
    for item in dbxref_str.split(","):
        item = item.strip()
        if ":" in item:
            k, v = item.split(":", 1)
            result[k.strip()] = v.strip()
    return result


def parse_refseq_gff(gff_paths, target_entrez_ids: Set[str]):
    """
    Parse one or more GFF/GFF3 files.
    Extracts mRNA → NM_xxx (transcript), and CDS → NP_xxx (protein) linked via Parent.
    
    Returns: {entrez_id -> set of (nm_versioned, np_versioned)}
    
    Strategy:
      - Pass 1: collect all mRNA rows → rna_id_to_nm: {rna_feature_id -> (nm_versioned, entrez_id)}
      - Pass 2: collect all CDS rows → parent_rna_id_to_np: {parent_rna_id -> np_versioned}
      - Combine: for each entrez_id, union (nm, np) pairs

    Handles both:
      .gff  format: Name=NM_xxx.v; Dbxref=GeneID:NNN,...
      .gff3 format: Name=NM_xxx.v; Dbxref=GeneID:NNN,... (similar)
    """
    # entrez_id -> set of (nm_versioned, np_versioned or "")
    result: Dict[str, Set[Tuple[str, str]]] = defaultdict(set)

    for gff_path in gff_paths:
        print(f"    Parsing: {os.path.basename(gff_path)} ...")
        # Pass 1: mRNA rows
        # feature_id (short like 'rna5' or 'rna-NM_xxx.v') -> nm_versioned
        feature_id_to_nm: Dict[str, str] = {}
        # nm_versioned -> entrez_id
        nm_to_entrez: Dict[str, str] = {}

        with open(gff_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                feature = parts[2]
                if feature != "mRNA":
                    continue
                attrs = _parse_gff_attrs(parts[8])
                name = attrs.get("Name", "")
                fid  = attrs.get("ID", "")
                dbxref_str = attrs.get("Dbxref", "")
                dbxref = _parse_dbxref(dbxref_str)
                entrez = dbxref.get("GeneID", "")

                # Only NM_ transcripts
                if not name.startswith("NM_"):
                    continue
                if entrez not in target_entrez_ids:
                    continue

                nm_versioned = name  # e.g. NM_007294.3
                # Store by the actual feature ID as used by CDS Parent=
                if fid:
                    feature_id_to_nm[fid] = nm_versioned
                # Also index by the rna-NM_xxx.v style (used in main GFF)
                feature_id_to_nm[f"rna-{nm_versioned}"] = nm_versioned
                nm_to_entrez[nm_versioned] = entrez

        # Pass 2: CDS rows to get NP_ linked to parent mRNA
        # Map: feature_id -> set of np_versioned
        parent_to_nps: Dict[str, Set[str]] = defaultdict(set)

        with open(gff_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] != "CDS":
                    continue
                attrs = _parse_gff_attrs(parts[8])
                parent = attrs.get("Parent", "")
                protein_id = attrs.get("protein_id", "")
                if not protein_id or not protein_id.startswith("NP_"):
                    continue
                if parent:
                    parent_to_nps[parent].add(protein_id)

        # Combine: for each known mRNA feature ID, merge NPs via parent lookup
        # Avoid duplicate (nm, np) from multiple feature_id aliases
        added: Set[Tuple[str, str, str]] = set()  # (entrez, nm, np)
        for fid, nm_versioned in feature_id_to_nm.items():
            entrez = nm_to_entrez.get(nm_versioned, "")
            if not entrez:
                continue
            nps = parent_to_nps.get(fid, set())
            if nps:
                for np in nps:
                    key = (entrez, nm_versioned, np)
                    if key not in added:
                        result[entrez].add((nm_versioned, np))
                        added.add(key)
            else:
                key = (entrez, nm_versioned, "")
                if key not in added:
                    result[entrez].add((nm_versioned, ""))
                    added.add(key)

    # Post-process: for each NM, if both (nm, "") and (nm, np) exist, drop the empty one
    for entrez in result:
        pairs = result[entrez]
        nms_with_np = {nm for nm, np in pairs if np}
        result[entrez] = {(nm, np) for nm, np in pairs if np or nm not in nms_with_np}

    return result


def parse_refseq_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Returns {np_versioned -> protein_sequence}
    Header example:
      >NP_570602.2 alpha-1-B glycoprotein precursor [Homo sapiens]
    """
    np_to_seq = {}
    current_id = None
    current_seq_parts = []

    _NP_RE = re.compile(r'(NP_\d+\.\d+)')
    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    np_to_seq[current_id] = "".join(current_seq_parts)
                # Handle both plain (>NP_xxx.v ...) and gi-format (>gi|N|ref|NP_xxx.v| ...)
                raw = line[1:].split()[0]
                m = _NP_RE.search(raw)
                current_id = m.group(1) if m else None
                current_seq_parts = []
            else:
                if current_id is not None:
                    current_seq_parts.append(line)
        if current_id is not None:
            np_to_seq[current_id] = "".join(current_seq_parts)

    return np_to_seq


def build_refseq_table(gene_set, gff_paths, fasta_paths, assembly: str):
    """
    Returns list of row dicts for the RefSeq transcript table.
    fasta_paths: list of protein FASTA paths; earlier paths take priority for sequence lookup.
    """
    # Build set of target entrez IDs
    entrez_to_gene = {}
    for gene in gene_set:
        eid = gene["entrez_id"]
        if eid:
            entrez_to_gene[eid] = gene

    target_entrez_ids = set(entrez_to_gene.keys())
    print(f"  Parsing {len(gff_paths)} GFF file(s) for {len(target_entrez_ids)} genes ...")

    entrez_to_transcripts = parse_refseq_gff(gff_paths, target_entrez_ids)

    print(f"  Parsing {len(fasta_paths)} RefSeq protein FASTA file(s) ...")
    np_to_seq: Dict[str, str] = {}
    for fasta_path in fasta_paths:
        print(f"    {os.path.basename(fasta_path)} ...")
        for np_id, seq in parse_refseq_fasta(fasta_path).items():
            if np_id not in np_to_seq:  # primary source has priority
                np_to_seq[np_id] = seq

    rows = []
    no_seq_count = 0

    for entrez, nm_np_pairs in sorted(entrez_to_transcripts.items()):
        gene = entrez_to_gene.get(entrez)
        if gene is None:
            continue
        for (nm_versioned, np_versioned) in sorted(nm_np_pairs):
            seq = np_to_seq.get(np_versioned, "") if np_versioned else ""
            if np_versioned and not seq:
                no_seq_count += 1
            rows.append({
                "gene_symbol":        gene["gene_symbol"],
                "entrez_id":          entrez,
                "hgnc_id":            gene["hgnc_id"],
                "ensembl_gene_id":    gene["ensembl_gene_id"],
                "refseq_transcript_id": nm_versioned,
                "refseq_protein_id":   np_versioned,
                "protein_sequence":    seq,
            })

    if no_seq_count:
        print(f"  ⚠️  {no_seq_count} NP IDs had no sequence in FASTA")

    return rows


def write_refseq_table(rows, path):
    cols = [
        "gene_symbol", "entrez_id", "hgnc_id", "ensembl_gene_id",
        "refseq_transcript_id", "refseq_protein_id", "protein_sequence"
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows)} rows → {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 60)
    print("Step 1: Collecting gene symbols from 4 sources ...")
    symbols = collect_gene_symbols()
    print(f"  Found {len(symbols)} unique gene symbols")

    print("\nStep 2: Loading HGNC lookup table ...")
    hgnc_path = fisoform("hgnc_complete_set_oct_07_2025.txt")
    symbol_to_info, alias_to_symbols, prev_to_symbols, hgnc_id_to_info = load_hgnc(hgnc_path)
    print(f"  Loaded {len(symbol_to_info)} HGNC entries ({len(hgnc_id_to_info)} with HGNC_ID)")

    print("\nStep 2b: Loading MANE fallback index ...")
    mane_path = fisoform("MANE.GRCh38.v1.2.summary.txt")
    mane_fallback = load_mane_fallback(mane_path, hgnc_id_to_info=hgnc_id_to_info)
    print(f"  Loaded {len(mane_fallback)} MANE symbols for fallback")

    print("\nStep 3: Building gene set ...")
    gene_set = build_gene_set(symbols, symbol_to_info, alias_to_symbols, prev_to_symbols,
                              mane_fallback=mane_fallback)
    write_gene_set(gene_set, outp("gene_set.tsv"))

    # Quick stats
    with_entrez = sum(1 for g in gene_set if g["entrez_id"])
    with_ensg   = sum(1 for g in gene_set if g["ensembl_gene_id"])
    print(f"  Genes with entrez_id: {with_entrez}/{len(gene_set)}")
    print(f"  Genes with ensembl_gene_id: {with_ensg}/{len(gene_set)}")

    # Build {current_sym -> [prev_sym, ...]} for the GTF fallback in build_ensembl_table
    sym_to_prev: Dict[str, list] = defaultdict(list)
    for prev_sym, current_syms in prev_to_symbols.items():
        for cur in current_syms:
            sym_to_prev[cur].append(prev_sym)

    print("\nStep 4: Building Ensembl GRCh37 transcript table ...")
    rows37 = build_ensembl_table(
        gene_set,
        finput("Homo_sapiens.GRCh37.87.gtf"),
        finput("Homo_sapiens.GRCh37.pep.all.fa"),
        "GRCh37",
        sym_to_prev=sym_to_prev,
    )
    write_ensembl_table(rows37, outp("ensembl_transcripts_grch37.tsv"))

    print("\nStep 5: Building Ensembl GRCh38 transcript table ...")
    rows38 = build_ensembl_table(
        gene_set,
        finput("Homo_sapiens.GRCh38.111.gtf"),
        finput("Homo_sapiens.GRCh38.pep.all.fa"),
        "GRCh38",
        sym_to_prev=sym_to_prev,
    )
    write_ensembl_table(rows38, outp("ensembl_transcripts_grch38.tsv"))

    print("\nStep 6: Building RefSeq GRCh37 transcript table ...")
    refseq37_rows = build_refseq_table(
        gene_set,
        gff_paths=[
            finput("GCF_000001405.25_GRCh37.p13_genomic.gff"),
            finput("ref_GRCh37.p10_top_level.gff3"),
            finput("ref_GRCh37.p13_top_level.gff3"),
        ],
        fasta_paths=[
            finput("GCF_000001405.25_GRCh37.p13_protein.faa"),
            finput("ref_GRCh37.p13_top_level.fa"),
            finput("ref_GRCh37.p10_top_level.fa"),
        ],
        assembly="GRCh37",
    )
    write_refseq_table(refseq37_rows, outp("refseq_transcripts_grch37.tsv"))

    print("\nStep 7: Building RefSeq GRCh38 transcript table ...")
    refseq38_rows = build_refseq_table(
        gene_set,
        gff_paths=[
            finput("GCF_000001405.40_GRCh38.p14_genomic.gff"),
            finput("GCF_000001405.38_GRCh38.p12_genomic.gff"),
        ],
        fasta_paths=[
            finput("GCF_000001405.40_GRCh38.p14_protein.faa"),
            finput("GCF_000001405.38_GRCh38.p12_protein.faa"),
        ],
        assembly="GRCh38",
    )
    write_refseq_table(refseq38_rows, outp("refseq_transcripts_grch38.tsv"))

    print("\n" + "=" * 60)
    print("Done. Output files in files/output/:")
    for fname in [
        "gene_set.tsv",
        "ensembl_transcripts_grch37.tsv",
        "ensembl_transcripts_grch38.tsv",
        "refseq_transcripts_grch37.tsv",
        "refseq_transcripts_grch38.tsv",
    ]:
        full = outp(fname)
        if os.path.exists(full):
            lines = sum(1 for _ in open(full)) - 1
            print(f"  {fname}: {lines} data rows")


if __name__ == "__main__":
    main()
