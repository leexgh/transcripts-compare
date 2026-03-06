import { useState, useEffect, useMemo } from 'react'
import type { GeneData, GeneEntry, FilterState } from '../lib/types'

export function useGeneData() {
  const [data, setData]     = useState<GeneData | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError]   = useState<string | null>(null)

  useEffect(() => {
    fetch('./data/gene_data.json')
      .then(r => {
        if (!r.ok) throw new Error(`HTTP ${r.status}`)
        return r.json()
      })
      .then((d: GeneData) => { setData(d); setLoading(false) })
      .catch(e => { setError(e.message); setLoading(false) })
  }, [])

  return { data, loading, error }
}

export function useFilteredGenes(genes: GeneEntry[], filters: FilterState): GeneEntry[] {
  return useMemo(() => {
    let result = genes
    const q = filters.search.trim().toLowerCase()
    if (q) {
      result = result.filter(g =>
        g.gene_symbol.toLowerCase().includes(q) ||
        g.best_match.grch37_enst.toLowerCase().includes(q) ||
        g.best_match.grch38_enst.toLowerCase().includes(q) ||
        g.best_match.grch37_nm.toLowerCase().includes(q) ||
        g.best_match.grch38_nm.toLowerCase().includes(q)
      )
    }
    if (filters.clinicalOnly)  result = result.filter(g => g.is_clinical)
    if (filters.germlineOnly)  result = result.filter(g => g.is_germline)
    if (filters.maneOnly)      result = result.filter(g => !!g.mane_grch38.enst)
    if (filters.oncokbOnly)    result = result.filter(g => g.is_oncokb)
    if (filters.mismatchOnly) {
      result = result.filter(g => {
        const s = g.similarities
        return (
          (s.grch37_enst_vs_grch37_refseq.pct !== null && s.grch37_enst_vs_grch37_refseq.pct < 100) ||
          (s.grch38_enst_vs_grch38_refseq.pct !== null && s.grch38_enst_vs_grch38_refseq.pct < 100) ||
          (s.grch37_enst_vs_grch38_enst.pct   !== null && s.grch37_enst_vs_grch38_enst.pct   < 100)
        )
      })
    }
    return result
  }, [genes, filters])
}
