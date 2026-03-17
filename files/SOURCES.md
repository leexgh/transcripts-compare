# Reference File Sources

This document lists where to download the large reference files used by the pipeline.

## Files

| File Name | Used In | Reference Genome | Version / Notes | Download URL |
|-----------|---------|-----------------|-----------------|--------------|
| `GCF_000001405.25_GRCh37.p13_protein.faa` | Find protein sequence | GRCh37 | Latest GRCh37 (p13) | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_protein.faa.gz) |
| `GCF_000001405.25_GRCh37.p13_genomic.gff` | Find NP and NM id | GRCh37 | Latest GRCh37 (p13) | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz) |
| `GCF_000001405.40_GRCh38.p14_protein.faa` | Find protein sequence | GRCh38 | 2025-08 (p14) | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz) |
| `GCF_000001405.40_GRCh38.p14_genomic.gff` | Find NP and NM id | GRCh38 | 2025-08 (p14) | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz) |
| `ref_GRCh37.p13_top_level.gff3` | Find NP and NM id | GRCh37 | Older annotation release 105 — needed to cover some RefSeq IDs not in p13 | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz) |
| `ref_GRCh37.p10_top_level.gff3` | Find NP and NM id | GRCh37 | Older annotation release 104 — needed to cover some RefSeq IDs | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.104/GFF/ref_GRCh37.p10_top_level.gff3.gz) |
| `GCF_000001405.38_GRCh38.p12_genomic.gff` | Find NP and NM id | GRCh38 | Older GRCh38 (p12) — only needed for 2 transcripts: `NM_203407.2` and `NM_001283009.1` | [Download](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gff.gz) |
| `hgnc_complete_set_oct_07_2025.txt` | HGNC gene table | — | Gene_Table_v7, HGNC Oct 07 2025 | [GitHub](https://github.com/cBioPortal/datahub-study-curation-tools/tree/master/gene-table-update/build-input-for-importer/gene-table-release-archives/Gene_Table_v7_HGNC_Oct_07_2025) |
| `Homo_sapiens.GRCh37.87.gtf` | Ensembl transcript IDs | GRCh37 | Ensembl release 87 (latest for GRCh37) | [Download](https://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz) |
| `Homo_sapiens.GRCh38.111.gtf` | Ensembl transcript IDs | GRCh38 | Ensembl release 111 | [Download](https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz) |
| `Homo_sapiens.GRCh38.pep.all.fa` | Ensembl protein sequences | GRCh38 | Ensembl release 111 | [Download](https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz) |
| `Homo_sapiens.GRCh37.pep.all.fa` | Ensembl protein sequences | GRCh37 | Ensembl release 111 via GRCh37 archive | [Download](https://ftp.ensembl.org/pub/grch37/release-111/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz) |

## Notes

- All files are distributed as `.gz` archives. Download and `gunzip` before use.
- `Homo_sapiens.GRCh38.116.gtf` (used in later pipeline runs) can be downloaded from:
  https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
- `Homo_sapiens.GRCh37.dna.primary_assembly.fa` and `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
  are the genome FASTA files. Download from Ensembl:
  - GRCh37: https://ftp.ensembl.org/pub/grch37/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
  - GRCh38: https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
