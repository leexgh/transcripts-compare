# Reference File Sources

This document lists where to download the large reference files used by the pipeline.
These files are **not** stored in the repository due to their size. Place them in `files/input/` before running the pipeline.

## RefSeq / NCBI Files (`files/input/`)

| File Name | Used In | Reference Genome | Version / Notes | Download |
|-----------|---------|-----------------|-----------------|----------|
| `GCF_000001405.25_GRCh37.p13_protein.faa` | Protein sequences | GRCh37 | Latest GRCh37 (p13) | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_protein.faa.gz) |
| `GCF_000001405.25_GRCh37.p13_genomic.gff` | NP and NM ID mapping | GRCh37 | Latest GRCh37 (p13) | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz) |
| `ref_GRCh37.p13_top_level.gff3` | NP and NM ID mapping | GRCh37 | Annotation release 105 — needed to cover some RefSeq IDs not in p13 | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz) |
| `ref_GRCh37.p10_top_level.gff3` | NP and NM ID mapping | GRCh37 | Annotation release 104 — needed to cover some additional RefSeq IDs | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.104/GFF/ref_GRCh37.p10_top_level.gff3.gz) |
| `protein.fa` (release 105) | Protein sequences | GRCh37 | Matches `ref_GRCh37.p13_top_level.gff3` — for older RefSeq transcripts | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/protein/protein.fa.gz) |
| `protein.fa` (release 104) | Protein sequences | GRCh37 | Matches `ref_GRCh37.p10_top_level.gff3` — for older RefSeq transcripts | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.104/protein/protein.fa.gz) |
| `GCF_000001405.40_GRCh38.p14_protein.faa` | Protein sequences | GRCh38 | 2025-08 (p14) | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz) |
| `GCF_000001405.40_GRCh38.p14_genomic.gff` | NP and NM ID mapping | GRCh38 | 2025-08 (p14) | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz) |
| `GCF_000001405.38_GRCh38.p12_genomic.gff` | NP and NM ID mapping | GRCh38 | Older GRCh38 (p12) — only needed for `NM_203407.2` and `NM_001283009.1` | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gff.gz) |
| `GCF_000001405.38_GRCh38.p12_protein.faa` | Protein sequences | GRCh38 | Older GRCh38 (p12) — supplementary protein sequences | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_protein.faa.gz) |

## Ensembl Files (`files/input/`)

| File Name | Used In | Reference Genome | Version / Notes | Download |
|-----------|---------|-----------------|-----------------|----------|
| `Homo_sapiens.GRCh38.87.gtf` | Ensembl transcript IDs | GRCh37 | Ensembl release 87 (latest for GRCh37 coordinate space) | [Download](https://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz) |
| `Homo_sapiens.GRCh37.pep.all.fa` | Ensembl protein sequences | GRCh37 | Ensembl release 111 via GRCh37 archive | [Download](https://ftp.ensembl.org/pub/grch37/release-111/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz) |
| `Homo_sapiens.GRCh38.111.gtf` | Ensembl transcript IDs | GRCh38 | Ensembl release 111 | [Download](https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz) |
| `Homo_sapiens.GRCh38.pep.all.fa` | Ensembl protein sequences | GRCh38 | Ensembl release 111 | [Download](https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz) |

## Isoform / Annotation Files (`files/isoform/`)

These smaller files are tracked in the repository under `files/isoform/`.

| File Name | Used In | Reference Genome | Version / Notes | Source |
|-----------|---------|-----------------|-----------------|--------|
| `isoform_overrides_at_mskcc_grch37.txt` | MSKCC GRCh37 isoform overrides | GRCh37 | — | [genome-nexus-importer](https://github.com/genome-nexus/genome-nexus-importer/blob/master/data/common_input/isoform_overrides_at_mskcc_grch37.txt) |
| `isoform_overrides_at_mskcc_grch38.txt` | MSKCC GRCh38 isoform overrides | GRCh38 | Draft | [leexgh/genome-nexus-importer](https://github.com/leexgh/genome-nexus-importer/blob/d8f8df3015d29803367c0251d94d9d945e0a1799/data/common_input/isoform_overrides_at_mskcc_grch38.txt) |
| `oncokb_grch37_canonical_transcripts.csv` / `oncokb_grch38_canonical_transcripts.csv` | OncoKB isoform overrides | Both | Provided via Slack — ideally sourced from [OncoKB API](https://www.oncokb.org/api/v1/utils/cancerGeneList.txt) | — |
| `MANE.GRCh38.v1.2.summary.txt` | MANE Select / Plus Clinical transcripts | GRCh38 | v1.2 (Ensembl 111) | [Download](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.summary.txt.gz) |
| `MANE_GRCh37_list_filtered.csv` | MANE GRCh37 liftover | GRCh37 | Filtered to Ensembl 111 transcripts only — sourced from [Ensembl Tark](https://tark.ensembl.org/web/mane_GRCh37_list/) | [Tark](https://tark.ensembl.org/web/mane_GRCh37_list/) |
| `hgnc_complete_set_oct_07_2025.txt` | HGNC gene table | — | Gene_Table_v7, HGNC Oct 07 2025 | [cBioPortal datahub](https://github.com/cBioPortal/datahub-study-curation-tools/tree/master/gene-table-update/build-input-for-importer/gene-table-release-archives/Gene_Table_v7_HGNC_Oct_07_2025) |
| `germline_panel_94.txt` | Germline gene panel | — | 94-gene panel | Internal |
| `Iv7_dmp_isoform_merged_overrides.txt` | Iv7 clinical isoform overrides | GRCh37 | DMP clinical panel | Internal |

## Notes

- All downloaded files are distributed as `.gz` archives. Download and `gunzip` before placing in `files/input/`.
- The two `protein.fa` files from release 104 and 105 should be renamed to avoid collision, e.g. `protein_r104.fa` and `protein_r105.fa`.
- OncoKB isoform files are currently shared via Slack; the plan is to migrate to the OncoKB API.
- The GRCh37 MANE liftover from Tark must be filtered to keep only transcripts matching Ensembl 111 versions.
