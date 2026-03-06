import React from 'react'
import type { SimilarityResult } from '../lib/types'

interface Props {
  sim:     SimilarityResult
  label?:  string
  onClick?: () => void
}

function pctClass(pct: number | null): string {
  if (pct === null) return 'bg-gray-100 text-gray-500'
  if (pct >= 100)   return 'bg-green-100 text-green-800'
  if (pct >= 95)    return 'bg-yellow-100 text-yellow-800'
  if (pct >= 90)    return 'bg-orange-100 text-orange-800'
  return 'bg-red-100 text-red-800'
}

export default function SimilarityBadge({ sim, label, onClick }: Props) {
  const { pct, diff_count } = sim
  const text = pct === null ? 'N/A' : `${pct.toFixed(pct === 100 ? 0 : 2)}%`

  return (
    <button
      type="button"
      onClick={onClick}
      title={label}
      className={`inline-flex items-center gap-1 px-2 py-0.5 rounded text-xs font-medium ${pctClass(pct)} ${onClick ? 'cursor-pointer hover:opacity-80 active:opacity-60' : 'cursor-default'}`}
    >
      {text}
      {diff_count > 0 && (
        <span className="opacity-70">({diff_count})</span>
      )}
    </button>
  )
}
