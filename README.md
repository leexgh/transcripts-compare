# Transcript Compare

A web app for comparing protein isoform transcripts across genome assemblies (GRCh37 vs GRCh38), data sources (Ensembl, RefSeq, MANE, OncoKB, MSKCC), and clinical panels.

**Live site:** [https://leexgh.github.io/transcripts-compare/](https://leexgh.github.io/transcripts-compare/)

## Features

### Transcript Table

The main view shows a filterable, sortable table of genes with their selected transcripts across four columns:

| Column | Description |
|--------|-------------|
| GRCh37 ENST | Ensembl transcript (GRCh37 assembly) |
| GRCh38 ENST | Ensembl transcript (GRCh38 assembly) |
| GRCh37 NM | RefSeq transcript (GRCh37 assembly) |
| GRCh38 NM | RefSeq transcript (GRCh38 assembly) |

Each transcript column has a **dropdown** — click to switch between all available transcripts for that gene. When you change a selection, similarity scores update in real time.

**Similarity columns** show pairwise protein sequence identity between selected transcripts (e.g., 37E/38E compares GRCh37 Ensembl vs GRCh38 Ensembl). Click any similarity badge to open the **diff viewer**.

**Badges:**
- `GL` = germline panel gene
- `CL` = clinical panel gene (Iv7)
- `Select` / `Plus Clinical` = MANE status

### Filters

- **Search** — filter by gene symbol or transcript ID
- **Assembly toggle** — show Both / GRCh37 only / GRCh38 only
- **Clinical only** — genes in the Iv7 clinical panel
- **Germline only** — genes in the germline panel
- **MANE only** — genes with a MANE Select or Plus Clinical transcript
- **Mismatches only** — genes where at least one similarity < 100%
- **OncoKB gene only** — genes in the OncoKB isoform list

### Diff Viewer

Click any similarity badge to open a side-by-side protein sequence alignment. Differences are color-coded:
- **Red** — substitutions (mismatched amino acids)
- **Blue** — insertions
- **Gray** — deletions (gaps)

### List Compare

Compare transcript IDs across different data sources side by side. Select two collections (e.g., OncoKB vs MANE) and see which genes match, differ, or are unique to each source. Includes similarity scores when both sources have sequences.

### Entrez Gene ID

Each gene row shows its NCBI Entrez Gene ID (when available), linked to the NCBI Gene page.

## Data Sources

### Isoform / Annotation Files (tracked in `files/isoform/`)

| Source | Description |
|--------|-------------|
| MSKCC isoform overrides (GRCh37) | Curated GRCh37 transcripts from [genome-nexus-importer](https://github.com/genome-nexus/genome-nexus-importer/blob/master/data/common_input/isoform_overrides_at_mskcc_grch37.txt) |
| MSKCC isoform overrides (GRCh38) | Draft GRCh38 overrides |
| Iv7 clinical overrides | DMP clinical panel RefSeq IDs (`Iv7_dmp_isoform_merged_overrides.txt`) |
| OncoKB isoform | GRCh37 + GRCh38 canonical transcripts (provided via Slack; see [OncoKB API](https://www.oncokb.org/api/v1/utils/cancerGeneList.txt)) |
| MANE GRCh38 | MANE Select / Plus Clinical — [v1.2](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.summary.txt.gz) (Ensembl 111) |
| MANE GRCh37 | GRCh37 liftover from [Ensembl Tark](https://tark.ensembl.org/web/mane_GRCh37_list/), filtered to Ensembl 111 transcripts |
| HGNC gene table | Gene_Table_v7, HGNC Oct 07 2025 via [cBioPortal datahub](https://github.com/cBioPortal/datahub-study-curation-tools/tree/master/gene-table-update/build-input-for-importer/gene-table-release-archives/Gene_Table_v7_HGNC_Oct_07_2025) |

### Large Reference Files (not tracked — see `files/SOURCES.md`)

| Source | Assembly | Description |
|--------|----------|-------------|
| RefSeq GRCh37 p13 GFF + FAA | GRCh37 | NM/NP ID mapping and protein sequences (latest GRCh37) |
| RefSeq annotation release 105 GFF3 + protein.fa | GRCh37 | Older release — covers additional RefSeq IDs |
| RefSeq annotation release 104 GFF3 + protein.fa | GRCh37 | Older release — covers additional RefSeq IDs |
| RefSeq GRCh38 p14 GFF + FAA | GRCh38 | NM/NP ID mapping and protein sequences (2025-08) |
| RefSeq GRCh38 p12 GFF + FAA | GRCh38 | Older release — only needed for `NM_203407.2` and `NM_001283009.1` |
| Ensembl release 87 GTF | GRCh37 | Ensembl transcript IDs for GRCh37 coordinate space |
| Ensembl release 111 GRCh37 pep FASTA | GRCh37 | Ensembl protein sequences |
| Ensembl release 111 GTF | GRCh38 | Ensembl transcript IDs |
| Ensembl release 111 GRCh38 pep FASTA | GRCh38 | Ensembl protein sequences |

See [`files/SOURCES.md`](files/SOURCES.md) for full download URLs.

## Development

### Prerequisites

- Node.js 18+
- npm

### Setup

```bash
git clone https://github.com/leexgh/transcripts-compare.git
cd transcripts-compare
npm install
```

### Run locally

```bash
npm run dev
```

Opens at [http://localhost:5173](http://localhost:5173).

### Build

```bash
npm run build
```

Output goes to `dist/`.

### Project structure

```
src/
  App.tsx                       # Router and main layout
  components/
    TranscriptTable.tsx          # Transcript Table tab
    FilterBar.tsx                # Search and filter controls
    TranscriptCompareTab.tsx     # Multiple Transcripts Compare tab
    ProteinComparePage.tsx       # Pairwise protein alignment page
    ListCompareTab.tsx           # List Compare tab
    DiffViewer.tsx               # Sequence diff viewer
    TranscriptDropdown.tsx       # Per-cell transcript selector
    SimilarityBadge.tsx          # Color-coded similarity display
    ComparePanel.tsx             # Multi-transcript comparison panel
    SequenceModal.tsx            # Full sequence display
  hooks/
    useGeneData.ts               # Data loading and filtering
    useCompareSelection.ts       # Comparison selection state
  lib/
    types.ts                     # TypeScript interfaces
    similarity.ts                # Client-side similarity computation
    alignment.ts                 # Needleman-Wunsch alignment
pipeline/
  build_gene_set.py             # Builds RefSeq/Ensembl TSVs from GFF + FASTA
  build_data.py                 # Builds gene_data.json from all sources
files/
  isoform/                      # Small annotation files (tracked)
  input/                        # Large reference files (not tracked, see SOURCES.md)
  output/                       # Generated TSVs (not tracked)
  SOURCES.md                    # Download URLs for large reference files
public/
  data/
    gene_data.json              # Pre-computed gene data (generated by pipeline)
```

## Deployment

The site auto-deploys to GitHub Pages on every push to `main` via GitHub Actions. No manual steps needed — just push your changes.

## Contributing

1. Fork this repo
2. Create a feature branch: `git checkout -b my-feature`
3. Make your changes and test locally with `npm run dev`
4. Commit and push to your fork
5. Open a Pull Request

Once merged to `main`, changes deploy automatically.
