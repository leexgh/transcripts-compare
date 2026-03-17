import React, { useMemo, useState } from 'react'
import { useGeneData } from '../hooks/useGeneData'
import { align, buildDiffBlocks } from '../lib/alignment'

// ---------------------------------------------------------------------------
// URL params: ?gene=PIK3CA&txs=ENST1,NM1&labels=MSKCC37,ClinRef
// ---------------------------------------------------------------------------
function parseParams() {
  const hash = window.location.hash  // e.g. "#/compare-protein?gene=...&txs=..."
  const qmark = hash.indexOf('?')
  if (qmark === -1) return { gene: '', txs: [], labels: [] }
  const raw = hash.slice(qmark + 1)
  const p = new URLSearchParams(raw)
  const txs = (p.get('txs') ?? '').split(',').map(s => s.trim()).filter(Boolean)
  const labelsRaw = (p.get('labels') ?? '').split(',').map(s => s.trim())
  const labels = txs.map((_, i) => labelsRaw[i] ?? txs[i])
  return { gene: p.get('gene') ?? '', txs, labels }
}

// ---------------------------------------------------------------------------
// Pct color helper
// ---------------------------------------------------------------------------
function pctBg(pct: number | null): string {
  if (pct === null) return 'bg-gray-100 text-gray-400'
  if (pct === 100)  return 'bg-green-100 text-green-800'
  if (pct >= 99)    return 'bg-yellow-50 text-yellow-800'
  if (pct >= 95)    return 'bg-orange-50 text-orange-700'
  return 'bg-red-50 text-red-700'
}

// ---------------------------------------------------------------------------
// Alignment viewer
// ---------------------------------------------------------------------------
const ROW_WIDTH = 60

function AlignmentView({ seq1, seq2, label1, label2 }: {
  seq1: string; seq2: string; label1: string; label2: string
}) {
  if (!seq1 || !seq2) {
    return <p className="text-gray-400 text-sm">Sequence not available for one or both transcripts.</p>
  }

  const result = useMemo(() => align(seq1, seq2), [seq1, seq2])
  const blocks = useMemo(() => buildDiffBlocks(result.seq1Aligned, result.seq2Aligned), [result])

  const subs = blocks.filter(b => b.type === 'sub').length
  const ins  = blocks.filter(b => b.type === 'ins').length
  const dels = blocks.filter(b => b.type === 'del').length

  // Build rows of ROW_WIDTH aligned columns
  const rows: Array<{ a1: string[]; a2: string[]; startPos1: number; startPos2: number }> = []
  const n = result.seq1Aligned.length
  let pos1 = 0, pos2 = 0

  for (let start = 0; start < n; start += ROW_WIDTH) {
    const chunk1 = result.seq1Aligned.slice(start, start + ROW_WIDTH).split('')
    const chunk2 = result.seq2Aligned.slice(start, start + ROW_WIDTH).split('')
    rows.push({ a1: chunk1, a2: chunk2, startPos1: pos1 + 1, startPos2: pos2 + 1 })
    for (let i = start; i < Math.min(start + ROW_WIDTH, n); i++) {
      if (result.seq1Aligned[i] !== '-') pos1++
      if (result.seq2Aligned[i] !== '-') pos2++
    }
  }

  function charColor(c1: string, c2: string, which: 1 | 2): string {
    const c = which === 1 ? c1 : c2
    if (c === '-') return 'text-red-500 bg-red-50'
    if (c1 === c2) return 'text-gray-700'
    return 'text-amber-700 bg-amber-50 font-bold'
  }

  return (
    <div>
      {/* Summary */}
      <div className="flex gap-4 mb-3 text-xs text-gray-600">
        <span className="font-medium">{result.pct.toFixed(2)}% identity</span>
        <span>{result.diffCount} positions differ</span>
        {subs > 0 && <span className="text-amber-600">{subs} substitution{subs !== 1 ? 's' : ''}</span>}
        {ins  > 0 && <span className="text-green-600">{ins} insertion{ins !== 1 ? 's' : ''}</span>}
        {dels > 0 && <span className="text-red-600">{dels} deletion{dels !== 1 ? 's' : ''}</span>}
      </div>

      {/* Alignment rows */}
      <div className="font-mono text-xs overflow-x-auto">
        {rows.map((row, ri) => (
          <div key={ri} className="mb-3">
            {/* Row 1 */}
            <div className="flex items-baseline gap-2">
              <span className="text-gray-400 w-32 shrink-0 truncate">{label1}</span>
              <span className="text-gray-400 w-8 text-right">{row.startPos1}</span>
              <span>
                {row.a1.map((c, i) => (
                  <span key={i} className={charColor(c, row.a2[i], 1)}>{c}</span>
                ))}
              </span>
            </div>
            {/* Row 2 */}
            <div className="flex items-baseline gap-2">
              <span className="text-gray-400 w-32 shrink-0 truncate">{label2}</span>
              <span className="text-gray-400 w-8 text-right">{row.startPos2}</span>
              <span>
                {row.a2.map((c, i) => (
                  <span key={i} className={charColor(row.a1[i], c, 2)}>{c}</span>
                ))}
              </span>
            </div>
          </div>
        ))}
      </div>

      {/* Diff table */}
      {blocks.length > 0 && (
        <details className="mt-4">
          <summary className="text-xs text-blue-600 cursor-pointer hover:underline">
            Show diff table ({blocks.length} events)
          </summary>
          <table className="mt-2 text-xs border-collapse w-full max-w-lg">
            <thead>
              <tr className="bg-gray-50">
                <th className="text-left px-2 py-1 border text-gray-600">Pos ({label1})</th>
                <th className="text-left px-2 py-1 border text-gray-600">Pos ({label2})</th>
                <th className="text-left px-2 py-1 border text-gray-600">{label1}</th>
                <th className="text-left px-2 py-1 border text-gray-600">{label2}</th>
                <th className="text-left px-2 py-1 border text-gray-600">Type</th>
              </tr>
            </thead>
            <tbody>
              {blocks.map((b, i) => (
                <tr key={i} className="border-b">
                  <td className="px-2 py-1 border">{b.pos1 + 1}</td>
                  <td className="px-2 py-1 border">{b.pos2 + 1}</td>
                  <td className="px-2 py-1 border font-mono">{b.seq1chars || '—'}</td>
                  <td className="px-2 py-1 border font-mono">{b.seq2chars || '—'}</td>
                  <td className="px-2 py-1 border capitalize">{b.type}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </details>
      )}
    </div>
  )
}

// ---------------------------------------------------------------------------
// Heatmap
// ---------------------------------------------------------------------------
function Heatmap({ txIds, labels, seqIndex, onSelect }: {
  txIds: string[]
  labels: string[]
  seqIndex: Map<string, string>
  onSelect: (i: number, j: number) => void
}) {
  const n = txIds.length
  const matrix = useMemo(() => {
    const m: Array<Array<number | null>> = []
    for (let i = 0; i < n; i++) {
      m.push([])
      for (let j = 0; j < n; j++) {
        if (i === j) { m[i].push(100); continue }
        const s1 = seqIndex.get(txIds[i]) ?? ''
        const s2 = seqIndex.get(txIds[j]) ?? ''
        if (!s1 || !s2) { m[i].push(null); continue }
        m[i].push(align(s1, s2).pct)
      }
    }
    return m
  }, [txIds, seqIndex])

  return (
    <div>
      <table className="border-collapse text-xs">
        <thead>
          <tr>
            <th />
            {labels.map((l, j) => (
              <th key={j} className="px-2 py-1 text-gray-600 font-medium max-w-24 truncate">{l}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {matrix.map((row, i) => (
            <tr key={i}>
              <td className="px-2 py-1 text-gray-600 font-medium text-right">{labels[i]}</td>
              {row.map((pct, j) => (
                <td
                  key={j}
                  className={`px-3 py-2 text-center cursor-pointer border ${pctBg(pct)} ${i === j ? 'opacity-50' : 'hover:ring-2 hover:ring-blue-300'}`}
                  onClick={() => i !== j && onSelect(i, j)}
                  title={`${labels[i]} vs ${labels[j]}: ${pct !== null ? pct.toFixed(2) + '%' : 'N/A'}`}
                >
                  {pct !== null ? (pct === 100 ? '100%' : `${pct.toFixed(2)}%`) : 'N/A'}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}

// ---------------------------------------------------------------------------
// Main page
// ---------------------------------------------------------------------------
export default function ProteinComparePage() {
  const { gene, txs, labels } = useMemo(parseParams, [])
  const { data, loading, error } = useGeneData()
  const [selected, setSelected] = useState<[number, number]>([0, Math.min(1, txs.length - 1)])

  const seqIndex = useMemo(() => {
    const m = new Map<string, string>()
    if (!data) return m
    for (const g of data.genes)
      for (const tx of g.transcripts)
        if (tx.sequence) m.set(tx.id, tx.sequence)
    return m
  }, [data])

  if (!txs.length) return (
    <div className="p-8 text-gray-500">No transcripts specified in URL.</div>
  )

  if (loading) return (
    <div className="p-8 text-gray-500">Loading sequence data…</div>
  )
  if (error) return (
    <div className="p-8 text-red-500">Error: {error}</div>
  )

  const [i, j] = selected
  const seq1 = seqIndex.get(txs[i]) ?? ''
  const seq2 = seqIndex.get(txs[j]) ?? ''

  return (
    <div className="p-6 max-w-5xl mx-auto">
      <h1 className="text-xl font-bold text-gray-800 mb-1">
        Protein Compare — {gene}
      </h1>
      <p className="text-xs text-gray-400 mb-6 font-mono break-all">
        {window.location.href}
      </p>

      {/* Transcript list */}
      <div className="mb-6">
        <h2 className="text-sm font-semibold text-gray-600 mb-2">Transcripts</h2>
        <div className="flex flex-wrap gap-2">
          {txs.map((tx, idx) => {
            const seq = seqIndex.get(tx)
            return (
              <span key={idx} className="inline-flex flex-col text-xs border rounded px-2 py-1 bg-white">
                <span className="font-medium text-blue-700">{labels[idx]}</span>
                <span className="font-mono text-gray-500">{tx}</span>
                <span className="text-gray-400">{seq ? `${seq.length} aa` : 'no sequence'}</span>
              </span>
            )
          })}
        </div>
      </div>

      {/* Heatmap */}
      <div className="mb-6">
        <h2 className="text-sm font-semibold text-gray-600 mb-2">Pairwise Similarity</h2>
        <p className="text-xs text-gray-400 mb-2">Click a cell to view alignment below.</p>
        <Heatmap
          txIds={txs}
          labels={labels}
          seqIndex={seqIndex}
          onSelect={(a, b) => setSelected([a, b])}
        />
      </div>

      {/* Alignment view */}
      <div>
        <h2 className="text-sm font-semibold text-gray-600 mb-2">
          Alignment: <span className="text-blue-600">{labels[i]}</span> vs <span className="text-blue-600">{labels[j]}</span>
        </h2>
        <div className="text-xs text-gray-400 mb-3 font-mono">
          <div>{txs[i]}</div>
          <div>{txs[j]}</div>
        </div>
        <AlignmentView
          seq1={seq1}
          seq2={seq2}
          label1={labels[i]}
          label2={labels[j]}
        />
      </div>
    </div>
  )
}
