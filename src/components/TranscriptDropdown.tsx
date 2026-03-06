import React, { useState, useRef, useEffect } from 'react'
import type { Transcript } from '../lib/types'

interface Props {
  current:     string
  transcripts: Transcript[]
  type:        'ensembl' | 'refseq'
  assembly:    'GRCh37' | 'GRCh38'
  onChange:    (t: Transcript) => void
}

export default function TranscriptDropdown({ current, transcripts, type, assembly, onChange }: Props) {
  const [open, setOpen] = useState(false)
  const ref = useRef<HTMLDivElement>(null)

  const alts    = transcripts.filter(t => t.type === type && t.assembly === assembly)
  const primary = alts.find(t => t.id === current)

  useEffect(() => {
    function handler(e: MouseEvent) {
      if (ref.current && !ref.current.contains(e.target as Node)) setOpen(false)
    }
    document.addEventListener('mousedown', handler)
    return () => document.removeEventListener('mousedown', handler)
  }, [])

  if (!current) {
    return <span className="text-gray-300 text-xs">—</span>
  }

  const tooltipText = primary?.source
    ? `Source: ${primary.source}${primary.is_primary ? '' : '\n(alternative)'}`
    : undefined

  if (alts.length <= 1) {
    return (
      <span className="font-mono text-xs" title={tooltipText}>
        {current}
      </span>
    )
  }

  return (
    <div className="relative inline-block" ref={ref}>
      <button
        onClick={() => setOpen(!open)}
        title={tooltipText}
        className="font-mono text-xs text-blue-600 hover:underline flex items-center gap-0.5 group"
      >
        {current}
        <span className="text-gray-400 group-hover:text-gray-600">▾</span>
      </button>
      {open && (
        <div className="absolute z-30 mt-1 left-0 bg-white border rounded shadow-lg min-w-max max-h-56 overflow-y-auto">
          {alts.map(t => (
            <button
              key={t.id}
              onClick={() => { onChange(t); setOpen(false) }}
              title={t.source ? `Source: ${t.source}` : undefined}
              className={`block w-full text-left px-3 py-1.5 text-xs hover:bg-gray-50 ${t.id === current ? 'bg-blue-50 text-blue-700' : ''}`}
            >
              <span className="font-mono">{t.id}</span>
              <span className="text-gray-400 ml-1">({t.length} AA)</span>
              {t.source && (
                <div className="text-[10px] text-gray-400 mt-0.5">{t.source}</div>
              )}
            </button>
          ))}
        </div>
      )}
    </div>
  )
}
