#!/usr/bin/env python3
"""
One-time script to fix outdated gene symbols in source data files.

Replacements (gene column only):
  WHSC1  -> NSD2   (Iv7 col 1)
  FTSJD1 -> CMTR2  (Iv7 col 1)
  STK19  -> WHR1   (Iv7 col 1, oncokb.csv col 2, oncokb_isoform.tsv col 0, MANE summary col 3)
  EIF2C1 -> AGO1   (Iv7 col 1)
  PAK7   -> PAK5   (Iv7 col 1)
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(ROOT, "files")

UPDATES = {
    "WHSC1":  "NSD2",
    "FTSJD1": "CMTR2",
    "NUT":    "NUTM1",
    "STK19":  "WHR1",
    "EIF2C1": "AGO1",
    "PAK7":   "PAK5",
}


def fix_tsv_column(filepath, col_idx, sep="\t"):
    """Replace old gene names in a specific column of a delimited file."""
    changed = 0
    lines = []
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            raw = line.rstrip("\n")
            parts = raw.split(sep)
            if col_idx < len(parts) and parts[col_idx].strip() in UPDATES:
                old = parts[col_idx].strip()
                # Preserve any surrounding whitespace in the field
                parts[col_idx] = parts[col_idx].replace(old, UPDATES[old])
                changed += 1
                lines.append(sep.join(parts) + "\n")
            else:
                lines.append(line)
    if changed:
        with open(filepath, "w", encoding="utf-8") as fh:
            fh.writelines(lines)
    return changed


def main():
    # 1. Iv7_dmp_isoform_merged_overrides.txt — gene is tab col 1
    iv7 = os.path.join(FILES, "Iv7_dmp_isoform_merged_overrides.txt")
    n = fix_tsv_column(iv7, col_idx=1, sep="\t")
    print(f"Iv7: {n} replacements")

    # 2. oncokb.csv — hugo_symbol is comma col 2
    oncokb_csv = os.path.join(FILES, "oncokb.csv")
    n = fix_tsv_column(oncokb_csv, col_idx=2, sep=",")
    print(f"oncokb.csv: {n} replacements")

    # 3. oncokb_isoform.tsv — Hugo Symbol is tab col 0
    oncokb_iso = os.path.join(FILES, "oncokb_isoform.tsv")
    n = fix_tsv_column(oncokb_iso, col_idx=0, sep="\t")
    print(f"oncokb_isoform.tsv: {n} replacements")

    # 4. MANE.GRCh38.v1.2.summary.txt — symbol is tab col 3
    mane = os.path.join(FILES, "MANE.GRCh38.v1.2.summary.txt")
    n = fix_tsv_column(mane, col_idx=3, sep="\t")
    print(f"MANE summary: {n} replacements")

    print("\nDone. Verify with: grep -E 'WHSC1|FTSJD1|EIF2C1|PAK7' in gene columns.")
    print("Note: STK19->WHR1 may still appear in comments/notes columns — that's expected.")


if __name__ == "__main__":
    main()
