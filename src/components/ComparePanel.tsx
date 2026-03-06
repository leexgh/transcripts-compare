import React from 'react'
import { useNavigate } from 'react-router-dom'
import type { Transcript } from '../lib/types'
import { computeSimilarity } from '../lib/similarity'
import SimilarityBadge from './SimilarityBadge'

interface Props {
  selected: Transcript[]
  onClear:  () => void
}

export default function ComparePanel({ selected, onClear }: Props) {
  const navigate = useNavigate()

  if (selected.length === 0) return null

  return (
    <div className="fixed bottom-0 left-0 right-0 bg-white border-t shadow-lg z-20 p-4">
      <div className="flex items-start justify-between max-w-7xl mx-auto">
        <div>
          <div className="flex items-center gap-2 mb-2">
            <span className="font-medium text-sm">{selected.length} transcript{selected.length > 1 ? 's' : ''} selected</span>
            <button onClick={onClear} className="text-xs text-gray-400 hover:text-gray-600">Clear</button>
          </div>

          <div className="flex flex-wrap gap-2">
            {selected.map(t => (
              <span key={t.id} className="px-2 py-1 bg-blue-50 rounded text-xs font-mono">
                {t.id} <span className="text-gray-400">({t.assembly})</span>
              </span>
            ))}
          </div>

          {selected.length >= 2 && (
            <div className="mt-3 grid gap-1">
              {selected.map((a, i) =>
                selected.slice(i + 1).map(b => {
                  const sim = computeSimilarity(a.sequence, b.sequence)
                  return (
                    <div key={`${a.id}-${b.id}`} className="flex items-center gap-2 text-xs">
                      <span className="font-mono">{a.id}</span>
                      <span className="text-gray-400">vs</span>
                      <span className="font-mono">{b.id}</span>
                      <SimilarityBadge
                        sim={sim}
                        onClick={() => navigate(`/diff?a=${encodeURIComponent(a.id)}&b=${encodeURIComponent(b.id)}`)}
                      />
                    </div>
                  )
                })
              )}
            </div>
          )}
        </div>

        {selected.length === 2 && (
          <button
            onClick={() => navigate(`/diff?a=${encodeURIComponent(selected[0].id)}&b=${encodeURIComponent(selected[1].id)}`)}
            className="px-4 py-2 bg-blue-600 text-white text-sm rounded hover:bg-blue-700 ml-4"
          >
            View Diff
          </button>
        )}
      </div>
    </div>
  )
}
