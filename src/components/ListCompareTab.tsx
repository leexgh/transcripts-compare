import React, { useState, useMemo, useCallback } from 'react'
import type { Collections, GeneData, SimilarityResult } from '../lib/types'
import { computeSimilarity } from '../lib/similarity'
import SimilarityBadge from './SimilarityBadge'

type AssemblyChoice = 'GRCh37' | 'GRCh38'
type IdType = 'ensembl' | 'refseq'
type CollectionKey = 'oncokb' | 'clinical' | 'mane' | 'mskcc_isoform'

const COLLECTIONS: { key: CollectionKey; label: string; assemblies: AssemblyChoice[] }[] = [
  { key: 'oncokb', label: 'OncoKB', assemblies: ['GRCh38'] },
  { key: 'clinical', label: 'Clinical (Iv7)', assemblies: ['GRCh37'] },
  { key: 'mane', label: 'MANE', assemblies: ['GRCh38'] },
  { key: 'mskcc_isoform', label: 'MSKCC Isoform', assemblies: ['GRCh37'] },
]

type RowStatus = 'Same' | 'Different' | 'Only in A' | 'Only in B'

interface CompareRow {
  gene: string
  idA: string
  idB: string
  status: RowStatus
  similarity: SimilarityResult
}

function extractId(entry: Record<string, string> | undefined, assembly: AssemblyChoice, idType: IdType): string {
  if (!entry) return ''
  if (idType === 'ensembl') {
    return (assembly === 'GRCh37' ? entry.grch37_enst : entry.grch38_enst) ?? ''
  }
  return (assembly === 'GRCh37' ? entry.grch37_nm : entry.grch38_nm) ?? ''
}

function statusClass(s: RowStatus): string {
  switch (s) {
    case 'Same': return 'bg-green-100 text-green-700'
    case 'Different': return 'bg-yellow-100 text-yellow-700'
    case 'Only in A': return 'bg-blue-100 text-blue-700'
    case 'Only in B': return 'bg-purple-100 text-purple-700'
  }
}

interface Props { data: GeneData }

export default function ListCompareTab({ data }: Props) {
  const [assembly, setAssembly] = useState<AssemblyChoice>('GRCh38')
  const [idType, setIdType] = useState<IdType>('ensembl')
  const [colA, setColA] = useState<CollectionKey>('oncokb')
  const [colB, setColB] = useState<CollectionKey>('mane')
  const [statusFilter, setStatusFilter] = useState<RowStatus | 'all'>('all')

  const availableCols = COLLECTIONS.filter(c => c.assemblies.includes(assembly))

  const rows = useMemo<CompareRow[]>(() => {
    const A = data.collections[colA] as Record<string, Record<string, string>>
    const B = data.collections[colB] as Record<string, Record<string, string>>
    const genes = new Set([...Object.keys(A), ...Object.keys(B)])

    // Build a map: transcript ID -> protein sequence from gene data
    const seqMap = new Map<string, string>()
    for (const g of data.genes) {
      for (const t of g.transcripts) {
        if (t.sequence && !seqMap.has(t.id)) {
          seqMap.set(t.id, t.sequence)
        }
      }
    }

    const out: CompareRow[] = []
    for (const gene of Array.from(genes).sort()) {
      const idA = extractId(A[gene], assembly, idType)
      const idB = extractId(B[gene], assembly, idType)
      let status: RowStatus
      if (!idA && !idB) continue
      if (!idA) status = 'Only in B'
      else if (!idB) status = 'Only in A'
      else if (idA === idB) status = 'Same'
      else status = 'Different'
      const seqA = seqMap.get(idA) ?? ''
      const seqB = seqMap.get(idB) ?? ''
      const similarity = computeSimilarity(seqA, seqB)
      out.push({ gene, idA, idB, status, similarity })
    }
    return out
  }, [data, colA, colB, assembly, idType])

  const filtered = useMemo(() =>
    statusFilter === 'all' ? rows : rows.filter(r => r.status === statusFilter),
    [rows, statusFilter]
  )

  const download = useCallback(() => {
    const header = 'Gene\tCollection A\tCollection B\tStatus\tSimilarity\n'
    const content = filtered.map(r => {
      const simStr = r.similarity.pct === null ? 'N/A' : `${r.similarity.pct.toFixed(2)}%`
      return `${r.gene}\t${r.idA}\t${r.idB}\t${r.status}\t${simStr}`
    }).join('\n')
    const blob = new Blob([header + content], { type: 'text/tab-separated-values' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url; a.download = 'transcript_comparison.tsv'; a.click()
    URL.revokeObjectURL(url)
  }, [filtered])

  const counts = useMemo(() => {
    const c = { Same: 0, Different: 0, 'Only in A': 0, 'Only in B': 0 }
    for (const r of rows) c[r.status]++
    return c
  }, [rows])

  return (
    <div className="flex flex-col flex-1 overflow-hidden p-4">
      {/* Controls */}
      <div className="flex flex-wrap gap-4 mb-4">
        <div>
          <label className="block text-xs font-medium text-gray-600 mb-1">Assembly</label>
          <div className="flex border rounded overflow-hidden text-sm">
            {(['GRCh37', 'GRCh38'] as AssemblyChoice[]).map(a => (
              <button key={a} onClick={() => setAssembly(a)}
                className={`px-3 py-1.5 ${assembly === a ? 'bg-blue-600 text-white' : 'hover:bg-gray-50'}`}>{a}</button>
            ))}
          </div>
        </div>

        <div>
          <label className="block text-xs font-medium text-gray-600 mb-1">ID Type</label>
          <div className="flex border rounded overflow-hidden text-sm">
            {([['ensembl', 'Ensembl (ENST)'], ['refseq', 'RefSeq (NM)']] as [IdType, string][]).map(([k, l]) => (
              <button key={k} onClick={() => setIdType(k)}
                className={`px-3 py-1.5 ${idType === k ? 'bg-blue-600 text-white' : 'hover:bg-gray-50'}`}>{l}</button>
            ))}
          </div>
        </div>

        <div>
          <label className="block text-xs font-medium text-gray-600 mb-1">Collection A</label>
          <select value={colA} onChange={e => setColA(e.target.value as CollectionKey)}
            className="border rounded px-2 py-1.5 text-sm">
            {availableCols.map(c => <option key={c.key} value={c.key}>{c.label}</option>)}
          </select>
        </div>

        <div>
          <label className="block text-xs font-medium text-gray-600 mb-1">Collection B</label>
          <select value={colB} onChange={e => setColB(e.target.value as CollectionKey)}
            className="border rounded px-2 py-1.5 text-sm">
            {availableCols.map(c => <option key={c.key} value={c.key}>{c.label}</option>)}
          </select>
        </div>

        <div className="self-end">
          <button onClick={download} className="px-3 py-1.5 border rounded text-sm hover:bg-gray-50">
            Download TSV
          </button>
        </div>
      </div>

      {/* Summary chips */}
      <div className="flex gap-2 mb-3 flex-wrap">
        <button onClick={() => setStatusFilter('all')}
          className={`px-2 py-0.5 rounded text-xs border ${statusFilter === 'all' ? 'bg-gray-700 text-white' : 'hover:bg-gray-50'}`}>
          All ({rows.length})
        </button>
        {(['Same', 'Different', 'Only in A', 'Only in B'] as RowStatus[]).map(s => (
          <button key={s} onClick={() => setStatusFilter(s)}
            className={`px-2 py-0.5 rounded text-xs border ${statusFilter === s ? 'border-gray-500 font-semibold' : 'hover:bg-gray-50'}`}>
            {s} ({counts[s]})
          </button>
        ))}
      </div>

      {/* Table */}
      <div className="flex-1 overflow-auto border rounded">
        <table className="min-w-full text-sm">
          <thead className="bg-gray-50 sticky top-0">
            <tr>
              <th className="text-left px-3 py-2 text-xs font-medium text-gray-600 border-b">Gene</th>
              <th className="text-left px-3 py-2 text-xs font-medium text-gray-600 border-b">Collection A</th>
              <th className="text-left px-3 py-2 text-xs font-medium text-gray-600 border-b">Collection B</th>
              <th className="text-left px-3 py-2 text-xs font-medium text-gray-600 border-b">Status</th>
              <th className="text-left px-3 py-2 text-xs font-medium text-gray-600 border-b">Similarity</th>
            </tr>
          </thead>
          <tbody>
            {filtered.map(row => (
              <tr key={row.gene} className="border-b hover:bg-gray-50">
                <td className="px-3 py-2 font-semibold">{row.gene}</td>
                <td className="px-3 py-2 font-mono text-xs">{row.idA || '—'}</td>
                <td className="px-3 py-2 font-mono text-xs">{row.idB || '—'}</td>
                <td className="px-3 py-2">
                  <span className={`px-2 py-0.5 rounded text-xs font-medium ${statusClass(row.status)}`}>{row.status}</span>
                </td>
                <td className="px-3 py-2">
                  <SimilarityBadge sim={row.similarity} label={`${row.idA || '—'} vs ${row.idB || '—'}`} />
                </td>
              </tr>
            ))}
          </tbody>
        </table>
        {filtered.length === 0 && (
          <div className="p-8 text-center text-gray-400 text-sm">No results.</div>
        )}
      </div>
    </div>
  )
}
