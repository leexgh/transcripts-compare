import React, { useState, useMemo, useCallback } from 'react'
import type { GeneEntry, GeneData } from '../lib/types'
import { computeSimilarity } from '../lib/similarity'

// ---------------------------------------------------------------------------
// Column definitions
// ---------------------------------------------------------------------------
type ColumnKey =
  | 'mskcc37_enst' | 'mskcc37_nm' | 'clinical_refseq'
  | 'oncokb37_enst' | 'oncokb37_nm' | 'mane37_enst' | 'mane37_nm'
  | 'mskcc38_enst' | 'mskcc38_nm'
  | 'oncokb38_enst' | 'oncokb38_nm' | 'mane38_enst' | 'mane38_nm'

type AssemblyFilter = 'GRCh37' | 'GRCh38' | 'Both'

interface ColDef {
  key: ColumnKey
  label: string
  assembly: 'GRCh37' | 'GRCh38' | 'Both'
  optional?: boolean
}

const ALL_COLS: ColDef[] = [
  { key: 'mskcc37_enst',    label: 'MSKCC37 ENST',      assembly: 'GRCh37' },
  { key: 'mskcc37_nm',      label: 'MSKCC37 RefSeq',    assembly: 'GRCh37' },
  { key: 'clinical_refseq', label: 'Clinical RefSeq',   assembly: 'GRCh37' },
  { key: 'oncokb37_enst',   label: 'OncoKB37 ENST',     assembly: 'GRCh37' },
  { key: 'oncokb37_nm',     label: 'OncoKB37 RefSeq',   assembly: 'GRCh37' },
  { key: 'mane37_enst',     label: 'MANE37 ENST',       assembly: 'GRCh37' },
  { key: 'mane37_nm',       label: 'MANE37 RefSeq',     assembly: 'GRCh37', optional: true },
  { key: 'mskcc38_enst',    label: 'MSKCC38 ENST',      assembly: 'GRCh38' },
  { key: 'mskcc38_nm',      label: 'MSKCC38 RefSeq',    assembly: 'GRCh38' },
  { key: 'oncokb38_enst',   label: 'OncoKB38 ENST',     assembly: 'GRCh38' },
  { key: 'oncokb38_nm',     label: 'OncoKB38 RefSeq',   assembly: 'GRCh38' },
  { key: 'mane38_enst',     label: 'MANE38 ENST',       assembly: 'GRCh38' },
  { key: 'mane38_nm',       label: 'MANE38 RefSeq',     assembly: 'GRCh38', optional: true },
]

const DEFAULT_COLS: Record<AssemblyFilter, ColumnKey[]> = {
  GRCh37: ['mskcc37_enst', 'clinical_refseq', 'oncokb37_nm', 'mane37_enst'],
  GRCh38: ['oncokb38_enst', 'oncokb38_nm', 'mane38_enst'],
  Both:   ['clinical_refseq', 'mskcc37_enst', 'mane38_enst', 'oncokb38_enst'],
}

// ---------------------------------------------------------------------------
// Value lookup
// ---------------------------------------------------------------------------
function getColValue(gene: GeneEntry, data: GeneData, key: ColumnKey): string {
  const c = data.collections
  switch (key) {
    case 'mskcc37_enst':    return c.mskcc_isoform[gene.gene_symbol]?.grch37_enst ?? ''
    case 'mskcc37_nm':      return c.mskcc_isoform[gene.gene_symbol]?.grch37_nm   ?? ''
    case 'clinical_refseq': return c.clinical[gene.gene_symbol]?.isoform_override  ?? c.clinical[gene.gene_symbol]?.grch37_nm ?? ''
    case 'oncokb37_enst':   return c.oncokb_isoform[gene.gene_symbol]?.grch37_enst ?? ''
    case 'oncokb37_nm':     return c.oncokb_isoform[gene.gene_symbol]?.grch37_nm   ?? ''
    case 'mane37_enst':     return gene.mane_grch37_enst                            ?? ''
    case 'mane37_nm':       return c.mane[gene.gene_symbol]?.grch38_nm             ?? ''
    case 'mskcc38_enst':    return c.mskcc_grch38_isoform[gene.gene_symbol]?.grch38_enst ?? ''
    case 'mskcc38_nm':      return c.mskcc_grch38_isoform[gene.gene_symbol]?.grch38_nm   ?? ''
    case 'oncokb38_enst':   return c.oncokb_isoform[gene.gene_symbol]?.grch38_enst ?? ''
    case 'oncokb38_nm':     return c.oncokb_isoform[gene.gene_symbol]?.grch38_nm   ?? ''
    case 'mane38_enst':     return c.mane[gene.gene_symbol]?.grch38_enst           ?? ''
    case 'mane38_nm':       return c.mane[gene.gene_symbol]?.grch38_nm             ?? ''
    default:                return ''
  }
}

function countResources(gene: GeneEntry, data: GeneData): number {
  return ALL_COLS.reduce((n, col) => n + (getColValue(gene, data, col.key) ? 1 : 0), 0)
}

// ---------------------------------------------------------------------------
// Sequence lookup (fuzzy: exact then base match)
// ---------------------------------------------------------------------------
function getSeqFuzzy(txId: string, gene: GeneEntry): string {
  if (!txId) return ''
  const exact = gene.transcripts.find(t => t.id === txId)
  if (exact?.sequence) return exact.sequence
  const base = txId.split('.')[0]
  return gene.transcripts.find(t => t.id.split('.')[0] === base)?.sequence ?? ''
}

function isAllIdenticalForCols(gene: GeneEntry, data: GeneData, cols: ColumnKey[]): boolean {
  const ids = cols.map(key => getColValue(gene, data, key)).filter(Boolean)
  if (ids.length < 2) return false
  const seqs = ids.map(id => getSeqFuzzy(id, gene)).filter(Boolean)
  if (seqs.length < 2) return false
  return seqs.every(s => s === seqs[0])
}

// ---------------------------------------------------------------------------
// CSV export
// ---------------------------------------------------------------------------
const BASE_URL = 'https://leexgh.github.io/transcripts-compare/'

function exportCSV(visibleGenes: GeneEntry[], data: GeneData, visibleCols: ColumnKey[]) {
  const colDefs = visibleCols.map(key => ALL_COLS.find(c => c.key === key)!)

  // All pairwise similarity column pairs
  const pairs: [number, number][] = []
  for (let i = 0; i < visibleCols.length; i++)
    for (let j = i + 1; j < visibleCols.length; j++)
      pairs.push([i, j])

  const headers = [
    'Gene', 'Clinical', 'Germline',
    ...colDefs.map(c => c.label),
    ...pairs.map(([i, j]) => `${colDefs[i].label} vs ${colDefs[j].label}`),
    'Link',
  ]

  const csvRows = visibleGenes.map(gene => {
    const vals  = visibleCols.map(key => getColValue(gene, data, key))
    const seqs  = vals.map(id => getSeqFuzzy(id, gene))

    const simCells = pairs.map(([i, j]) => {
      if (!vals[i] || !vals[j]) return ''
      const sim = computeSimilarity(seqs[i], seqs[j])
      return sim.pct === null ? 'N/A' : `${sim.pct.toFixed(sim.pct === 100 ? 0 : 2)}%`
    })

    const validIdxs  = vals.map((v, i) => v ? i : -1).filter(i => i >= 0)
    const validIds   = validIdxs.map(i => vals[i])
    const validLabels = validIdxs.map(i => colDefs[i].label)
    const link = validIds.length >= 2
      ? `${BASE_URL}#/compare-protein?gene=${encodeURIComponent(gene.gene_symbol)}&txs=${encodeURIComponent(validIds.join(','))}&labels=${encodeURIComponent(validLabels.join(','))}`
      : ''

    return [
      gene.gene_symbol,
      gene.is_clinical ? 'Yes' : 'No',
      gene.is_germline ? 'Yes' : 'No',
      ...vals,
      ...simCells,
      link,
    ]
  })

  const escape = (s: string) => `"${String(s ?? '').replace(/"/g, '""')}"`
  const csv = [headers, ...csvRows].map(r => r.map(escape).join(',')).join('\n')
  const blob = new Blob([csv], { type: 'text/csv' })
  const url  = URL.createObjectURL(blob)
  const a    = document.createElement('a')
  a.href = url
  a.download = 'transcript_compare.csv'
  a.click()
  URL.revokeObjectURL(url)
}

// ---------------------------------------------------------------------------
// Multi-compare cell
// ---------------------------------------------------------------------------
interface Entry { id: string; label: string }

interface MultiCompareCellProps {
  gene: GeneEntry
  data: GeneData
  visibleCols: ColumnKey[]
}

function MultiCompareCell({ gene, data, visibleCols }: MultiCompareCellProps) {
  const [entries, setEntries] = useState<Entry[]>(() =>
    visibleCols
      .map(key => ({
        id:    getColValue(gene, data, key),
        label: ALL_COLS.find(c => c.key === key)?.label ?? key,
      }))
      .filter(e => Boolean(e.id))
  )

  const allTxIds = useMemo(() =>
    Array.from(new Set(gene.transcripts.map(t => t.id))).sort()
  , [gene])

  const validEntries = entries.filter(e => Boolean(e.id))
  const validIds     = validEntries.map(e => e.id)

  const allIdentical = useMemo(() => {
    if (validIds.length < 2) return false
    const seqs = validIds.map(id => getSeqFuzzy(id, gene)).filter(Boolean)
    if (seqs.length < 2) return false
    return seqs.every(s => s === seqs[0])
  }, [validIds, gene])

  const handleCompare = useCallback(() => {
    if (validIds.length < 2) return
    const labels = validEntries.map(e => e.label || e.id)
    const url = `#/compare-protein?gene=${encodeURIComponent(gene.gene_symbol)}&txs=${encodeURIComponent(validIds.join(','))}&labels=${encodeURIComponent(labels.join(','))}`
    window.open(window.location.href.split('#')[0] + url, '_blank')
  }, [validEntries, gene])

  return (
    <div className="flex flex-col gap-1 min-w-[300px]">
      {entries.map((entry, idx) => (
        <div key={idx} className="flex items-center gap-0.5">
          <select
            value={entry.id}
            onChange={e => setEntries(prev => prev.map((v, i) =>
              i === idx ? { id: e.target.value, label: e.target.value } : v
            ))}
            className="text-xs border rounded px-1 py-0.5 font-mono w-[260px] truncate bg-white"
          >
            <option value="">—</option>
            {allTxIds.map(t => <option key={t} value={t}>{t}</option>)}
          </select>
          {entries.length > 2 && (
            <button
              onClick={() => setEntries(prev => prev.filter((_, i) => i !== idx))}
              className="text-gray-300 hover:text-red-400 text-xs leading-none"
              title="Remove"
            >×</button>
          )}
        </div>
      ))}
      <div className="flex items-center gap-1.5 mt-0.5">
        <button onClick={() => setEntries(prev => [...prev, { id: '', label: '' }])} className="text-xs text-blue-500 hover:underline">+ Add</button>
        <button
          onClick={handleCompare}
          disabled={validIds.length < 2}
          className="px-2 py-0.5 text-xs bg-blue-50 text-blue-600 border border-blue-200 rounded hover:bg-blue-100 disabled:opacity-40"
          title="Open protein compare"
        >⊞</button>
        {allIdentical && <span className="text-xs font-medium text-green-600">100%</span>}
      </div>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Main component
// ---------------------------------------------------------------------------
const PAGE_SIZE = 500

interface Props { data: GeneData }

export default function TranscriptCompareTab({ data }: Props) {
  const [assembly, setAssembly]       = useState<AssemblyFilter>('GRCh37')
  const [search, setSearch]           = useState('')
  const [visibleCols, setVisibleCols] = useState<ColumnKey[]>(DEFAULT_COLS['GRCh37'])
  const [maneLimit, setManeLimit]     = useState(0)
  const [mismatchOnly, setMismatchOnly] = useState(false)

  const handleAssembly = useCallback((a: AssemblyFilter) => {
    setAssembly(a)
    setVisibleCols(DEFAULT_COLS[a])
    setManeLimit(0)
  }, [])

  const availableCols = useMemo(() =>
    assembly === 'Both' ? ALL_COLS : ALL_COLS.filter(c => c.assembly === assembly)
  , [assembly])

  const addCol    = useCallback((key: ColumnKey) => {
    if (!visibleCols.includes(key)) setVisibleCols(prev => [...prev, key])
  }, [visibleCols])
  const removeCol = useCallback((key: ColumnKey) => setVisibleCols(prev => prev.filter(k => k !== key)), [])
  const swapCol   = useCallback((idx: number, newKey: ColumnKey) => setVisibleCols(prev => prev.map((k, i) => i === idx ? newKey : k)), [])

  const optionalCols = useMemo(() => availableCols.filter(c => !visibleCols.includes(c.key)), [availableCols, visibleCols])

  // ---------------------------------------------------------------------------
  // Gene lists
  // ---------------------------------------------------------------------------
  const q        = search.trim().toLowerCase()
  const allGenes = data.genes

  const filteredGenes = useMemo(() => {
    if (!q) return allGenes
    return allGenes.filter(g =>
      g.gene_symbol.toLowerCase().includes(q) ||
      (g.hgnc_id && g.hgnc_id.toLowerCase().includes(q)) ||
      (g.entrez_gene_id && g.entrez_gene_id.includes(q))
    )
  }, [allGenes, q])

  const curatedGenes = useMemo(() =>
    filteredGenes
      .filter(g => !g.mane_only)
      .sort((a, b) => countResources(b, data) - countResources(a, data) || a.gene_symbol.localeCompare(b.gene_symbol))
  , [filteredGenes, data])

  const maneOnlyGenes = useMemo(() => filteredGenes.filter(g => g.mane_only), [filteredGenes])

  const baseGenes = useMemo(() => {
    if (q) return filteredGenes
    return [...curatedGenes, ...maneOnlyGenes.slice(0, maneLimit)]
  }, [q, filteredGenes, curatedGenes, maneOnlyGenes, maneLimit])

  const visibleGenes = useMemo(() => {
    if (!mismatchOnly) return baseGenes
    return baseGenes.filter(g => !isAllIdenticalForCols(g, data, visibleCols))
  }, [baseGenes, mismatchOnly, data, visibleCols])

  const remainingMane = q ? 0 : Math.max(0, maneOnlyGenes.length - maneLimit)
  const sec2ColSpan   = visibleCols.length + (optionalCols.length > 0 ? 1 : 0)

  return (
    <div className="flex flex-col flex-1 overflow-hidden">
      {/* Top bar */}
      <div className="flex items-center gap-4 px-4 py-2 bg-white border-b">
        <div className="flex rounded border overflow-hidden text-sm">
          {(['GRCh37', 'GRCh38', 'Both'] as AssemblyFilter[]).map(a => (
            <button
              key={a}
              onClick={() => handleAssembly(a)}
              className={`px-3 py-1 ${assembly === a ? 'bg-blue-600 text-white' : 'bg-white text-gray-600 hover:bg-gray-50'}`}
            >
              {a}
            </button>
          ))}
        </div>

        <input
          type="text"
          value={search}
          onChange={e => { setSearch(e.target.value); setManeLimit(0) }}
          placeholder="Filter genes…"
          className="border rounded px-2 py-1 text-sm w-48 focus:outline-none focus:ring-1 focus:ring-blue-400"
        />

        <button
          onClick={() => setMismatchOnly(v => !v)}
          className={`px-3 py-1 text-sm rounded border ${mismatchOnly ? 'bg-orange-500 text-white border-orange-500' : 'bg-white text-gray-600 border-gray-300 hover:bg-gray-50'}`}
        >
          Mismatch only
        </button>

        <span className="text-xs text-gray-400 ml-auto">
          {visibleGenes.length} gene{visibleGenes.length !== 1 ? 's' : ''}
          {!q && remainingMane > 0 && ` (${remainingMane} MANE-only hidden)`}
        </span>

        <button
          onClick={() => exportCSV(visibleGenes, data, visibleCols)}
          className="px-3 py-1 text-sm rounded border bg-white text-gray-600 border-gray-300 hover:bg-gray-50"
          title="Download current table as CSV"
        >
          ↓ Export CSV
        </button>
      </div>

      {/* Table */}
      <div className="flex-1 overflow-auto">
        <table className="min-w-full text-sm border-collapse">
          <thead className="bg-gray-50 sticky top-0 z-10">
            <tr>
              <th className="px-3 py-1 text-xs font-semibold text-gray-500 border-b border-r bg-gray-100 text-center">Gene</th>
              <th colSpan={sec2ColSpan} className="px-3 py-1 text-xs font-semibold text-gray-500 border-b border-r bg-gray-100 text-center">
                Transcript IDs
              </th>
              <th className="px-3 py-1 text-xs font-semibold text-gray-500 border-b bg-gray-100 text-center">
                Compare
              </th>
            </tr>
            <tr>
              <th className="text-left px-3 py-2 text-xs font-medium text-gray-600 border-b border-r w-28 whitespace-nowrap">Gene</th>

              {visibleCols.map((key, idx) => (
                <th key={key} className="text-left px-2 py-1 text-xs font-medium text-gray-600 border-b">
                  <div className="flex items-center gap-1">
                    <select
                      value={key}
                      onChange={e => swapCol(idx, e.target.value as ColumnKey)}
                      className="text-xs border-0 bg-transparent font-medium text-gray-700 focus:outline-none cursor-pointer pr-1"
                    >
                      {availableCols.map(c => <option key={c.key} value={c.key}>{c.label}</option>)}
                    </select>
                    <button onClick={() => removeCol(key)} className="text-gray-300 hover:text-red-400 text-xs leading-none" title="Remove column">×</button>
                  </div>
                </th>
              ))}
              {optionalCols.length > 0 ? (
                <th className="px-2 py-1 border-b border-r">
                  <select
                    value=""
                    onChange={e => { if (e.target.value) addCol(e.target.value as ColumnKey) }}
                    className="text-xs border rounded px-1 py-0.5 text-gray-500 bg-white focus:outline-none"
                  >
                    <option value="">+ Add col</option>
                    {optionalCols.map(c => <option key={c.key} value={c.key}>{c.label}</option>)}
                  </select>
                </th>
              ) : (
                <th className="border-b border-r w-0 p-0" />
              )}

              <th className="text-left px-2 py-2 text-xs font-medium text-gray-600 border-b whitespace-nowrap min-w-[300px]">Multi Compare</th>
            </tr>
          </thead>

          <tbody>
            {visibleGenes.map(gene => {
              const values = visibleCols.map(key => getColValue(gene, data, key))
              return (
                <tr key={gene.gene_symbol} className="border-b hover:bg-gray-50">
                  <td className="px-3 py-2 align-top border-r">
                    <div className="font-semibold text-sm">{gene.gene_symbol}</div>
                    <div className="flex gap-1 mt-0.5">
                      {gene.is_germline && <span className="px-1 py-0.5 bg-purple-100 text-purple-700 text-[10px] rounded">Germline</span>}
                      {gene.is_clinical  && <span className="px-1 py-0.5 bg-blue-100 text-blue-700 text-[10px] rounded">Clinical</span>}
                      {gene.mane_only    && <span className="px-1 py-0.5 bg-gray-100 text-gray-500 text-[10px] rounded">MANE</span>}
                    </div>
                  </td>

                  {values.map((val, idx) => (
                    <td key={visibleCols[idx]} className="px-2 py-2 align-top font-mono text-xs">
                      {val ? <span className="text-gray-700">{val}</span> : <span className="text-gray-300">—</span>}
                    </td>
                  ))}
                  <td className="border-r w-0 p-0" />

                  <td className="px-2 py-2 align-top">
                    <MultiCompareCell key={visibleCols.join(',')} gene={gene} data={data} visibleCols={visibleCols} />
                  </td>
                </tr>
              )
            })}
          </tbody>
        </table>
      </div>

      {!q && remainingMane > 0 && (
        <div className="border-t py-2 text-center">
          <button
            onClick={() => setManeLimit(prev => prev + PAGE_SIZE)}
            className="text-sm text-blue-600 hover:underline"
          >
            Show more
          </button>
        </div>
      )}
    </div>
  )
}
