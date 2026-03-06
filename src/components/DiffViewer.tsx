import React, { useEffect, useState, useMemo } from 'react'
import { useSearchParams } from 'react-router-dom'
import type { GeneData, Transcript } from '../lib/types'
import { align, buildDiffBlocks } from '../lib/alignment'
import type { DiffBlock } from '../lib/alignment'

const CHARS_PER_ROW = 60

interface Props {
  data:   GeneData | null
  seq1?:  string
  seq2?:  string
  name1?: string
  name2?: string
}

function lookupSeq(data: GeneData | null, id: string): string {
  if (!data) return ''
  for (const gene of data.genes) {
    const t = gene.transcripts.find(t => t.id === id)
    if (t) return t.sequence
  }
  return ''
}

function charClass(c1: string, c2: string): string {
  if (c1 === '-' && c2 !== '-') return 'bg-blue-200 text-blue-900'
  if (c1 !== '-' && c2 === '-') return 'bg-gray-200 text-gray-600'
  if (c1 !== c2)                return 'bg-red-200 text-red-900 font-bold'
  return ''
}

function AlignedRow({ a1, a2, offset }: { a1: string; a2: string; offset: number }) {
  const pos = offset + 1
  return (
    <div className="font-mono text-xs mb-2">
      <div className="text-gray-400 mb-0.5 select-none">
        {String(pos).padStart(6, ' ')}
      </div>
      <div className="flex">
        <span className="text-gray-500 mr-2 select-none">Seq1:</span>
        {Array.from(a1).map((c, i) => (
          <span key={i} className={charClass(c, a2[i])}>{c}</span>
        ))}
      </div>
      <div className="flex">
        <span className="text-gray-400 mr-2 select-none">     </span>
        {Array.from(a1).map((c, i) => (
          <span key={i} className="text-gray-400">{c === a2[i] ? '|' : ' '}</span>
        ))}
      </div>
      <div className="flex">
        <span className="text-gray-500 mr-2 select-none">Seq2:</span>
        {Array.from(a2).map((c, i) => (
          <span key={i} className={charClass(a2[i], a1[i])}>{c}</span>
        ))}
      </div>
    </div>
  )
}

export default function DiffViewer({ data, seq1: propSeq1, seq2: propSeq2, name1: propName1, name2: propName2 }: Props) {
  const [searchParams] = useSearchParams()
  const id1   = searchParams.get('a') ?? ''
  const id2   = searchParams.get('b') ?? ''

  const seq1  = propSeq1 ?? lookupSeq(data, id1)
  const seq2  = propSeq2 ?? lookupSeq(data, id2)
  const name1 = propName1 ?? id1
  const name2 = propName2 ?? id2

  const result = useMemo(() => {
    if (!seq1 || !seq2) return null
    return align(seq1, seq2)
  }, [seq1, seq2])

  const blocks = useMemo(() => {
    if (!result) return []
    return buildDiffBlocks(result.seq1Aligned, result.seq2Aligned)
  }, [result])

  const rows = useMemo(() => {
    if (!result) return []
    const out: { a1: string; a2: string; offset: number }[] = []
    const n = result.seq1Aligned.length
    let pos1 = 0
    for (let i = 0; i < n; i += CHARS_PER_ROW) {
      const chunk1 = result.seq1Aligned.slice(i, i + CHARS_PER_ROW)
      const chunk2 = result.seq2Aligned.slice(i, i + CHARS_PER_ROW)
      out.push({ a1: chunk1, a2: chunk2, offset: pos1 })
      pos1 += chunk1.replace(/-/g, '').length
    }
    return out
  }, [result])

  if (!seq1 || !seq2) {
    return (
      <div className="p-8 text-center text-gray-500">
        {!seq1 && !seq2 ? 'No sequences provided.' : `Missing sequence for ${!seq1 ? name1 : name2}`}
      </div>
    )
  }

  return (
    <div className="p-4">
      <div className="mb-4">
        <h2 className="text-lg font-semibold">
          {name1} <span className="text-gray-400">vs</span> {name2}
          {result && (
            <span className="ml-3 text-base font-normal text-gray-600">
              — {result.pct.toFixed(2)}% identical ({result.diffCount} difference{result.diffCount !== 1 ? 's' : ''})
            </span>
          )}
        </h2>
        <div className="text-sm text-gray-500 mt-1">
          Seq1: {seq1.length} AA &nbsp;|&nbsp; Seq2: {seq2.length} AA
        </div>
      </div>

      {/* Legend */}
      <div className="flex gap-3 mb-4 text-xs">
        <span className="px-2 py-0.5 bg-red-200 rounded">Substitution</span>
        <span className="px-2 py-0.5 bg-blue-200 rounded">Insertion</span>
        <span className="px-2 py-0.5 bg-gray-200 rounded">Deletion</span>
      </div>

      {/* Aligned view */}
      <div className="bg-gray-50 rounded-lg p-4 overflow-x-auto mb-6">
        {rows.map((r, i) => (
          <AlignedRow key={i} a1={r.a1} a2={r.a2} offset={r.offset} />
        ))}
      </div>

      {/* Diff summary */}
      {blocks.length > 0 && (
        <div className="mt-4">
          <h3 className="font-medium mb-2 text-sm">Differences</h3>
          <div className="divide-y text-sm">
            {blocks.map((b, i) => (
              <div key={i} className="py-1.5 flex items-start gap-3">
                <span className={`px-2 py-0.5 rounded text-xs font-medium ${
                  b.type === 'sub' ? 'bg-red-100 text-red-700' :
                  b.type === 'ins' ? 'bg-blue-100 text-blue-700' :
                  'bg-gray-100 text-gray-600'
                }`}>{b.type}</span>
                <span>
                  {b.type === 'sub' && (
                    <><code className="font-mono">{b.seq1chars}</code> → <code className="font-mono">{b.seq2chars}</code> at position {b.pos1 + 1}</>
                  )}
                  {b.type === 'ins' && (
                    <><code className="font-mono">{b.seq2chars}</code> insertion at position {b.pos2 + 1} in Seq2</>
                  )}
                  {b.type === 'del' && (
                    <><code className="font-mono">{b.seq1chars}</code> deletion at position {b.pos1 + 1} in Seq1</>
                  )}
                </span>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  )
}
