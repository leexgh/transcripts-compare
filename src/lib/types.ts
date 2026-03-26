export interface SimilarityResult {
  pct:        number | null
  diff_count: number
}

export interface Similarities {
  grch37_enst_vs_grch37_refseq: SimilarityResult
  grch38_enst_vs_grch38_refseq: SimilarityResult
  grch37_enst_vs_grch38_enst:   SimilarityResult
  grch37_enst_vs_mane:          SimilarityResult
  grch38_enst_vs_mane:          SimilarityResult
  grch37_nm_vs_mane?:           SimilarityResult
}

export interface Transcript {
  id:         string
  type:       'ensembl' | 'refseq'
  assembly:   'GRCh37' | 'GRCh38'
  sequence:   string
  length:     number
  source:     string
  is_primary: boolean
  ccds?:      string
}

export interface BestMatch {
  grch37_enst: string
  grch38_enst: string
  grch37_nm:   string
  grch38_nm:   string
}

export interface ManeGRCh38 {
  enst:   string
  nm:     string
  status: string
}

export interface GeneEntry {
  gene_symbol:        string
  gene_id:            string
  hgnc_id:            string
  entrez_gene_id:     string
  is_germline:        boolean
  is_clinical:        boolean
  is_oncokb:          boolean
  mane_only:          boolean
  clinical_refseq:    string
  curated_note:       string
  mskcc38_enst:       string
  mskcc38_nm:         string
  mane_grch38:        ManeGRCh38
  mane_grch37_enst:   string
  best_match:         BestMatch
  oncokb_best_match?: BestMatch
  transcripts:        Transcript[]
  similarities:       Similarities
}

export interface CollectionEntry {
  grch38_enst?: string
  grch38_nm?:   string
  grch37_enst?: string
  grch37_nm?:   string
  status?:      string
}

export interface Collections {
  oncokb:                Record<string, CollectionEntry>
  oncokb_isoform:        Record<string, CollectionEntry>
  clinical:              Record<string, { grch37_nm: string; isoform_override?: string }>
  mane:                  Record<string, CollectionEntry>
  mskcc_isoform:         Record<string, CollectionEntry>
  mskcc_grch38_isoform:  Record<string, CollectionEntry>
}

export interface Metadata {
  generated_at:    string
  ensembl_version: number
  mane_version:    string
  gene_count:      number
}

export interface GeneData {
  metadata:    Metadata
  genes:       GeneEntry[]
  collections: Collections
}

export type AssemblyFilter = 'both' | 'GRCh37' | 'GRCh38'

export interface FilterState {
  assembly:      AssemblyFilter
  search:        string
  clinicalOnly:  boolean
  germlineOnly:  boolean
  maneOnly:      boolean
  mismatchOnly:  boolean
}
