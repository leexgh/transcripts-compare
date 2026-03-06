import React, { useCallback } from 'react'
import type { Transcript } from '../lib/types'

interface Props {
  transcript: Transcript | null
  onClose:    () => void
}

function formatFasta(id: string, seq: string): string {
  const lines = [`>${id}`]
  for (let i = 0; i < seq.length; i += 60) {
    lines.push(seq.slice(i, i + 60))
  }
  return lines.join('\n')
}

export default function SequenceModal({ transcript, onClose }: Props) {
  if (!transcript) return null

  const fasta = formatFasta(transcript.id, transcript.sequence)

  const copy = useCallback(() => {
    navigator.clipboard.writeText(fasta)
  }, [fasta])

  const download = useCallback(() => {
    const blob = new Blob([fasta], { type: 'text/plain' })
    const url  = URL.createObjectURL(blob)
    const a    = document.createElement('a')
    a.href     = url
    a.download = `${transcript.id}.fasta`
    a.click()
    URL.revokeObjectURL(url)
  }, [fasta, transcript.id])

  return (
    <div className="fixed inset-0 bg-black/40 z-50 flex items-center justify-center p-4" onClick={onClose}>
      <div className="bg-white rounded-xl shadow-2xl w-full max-w-2xl max-h-[80vh] flex flex-col" onClick={e => e.stopPropagation()}>
        <div className="flex items-center justify-between p-4 border-b">
          <div>
            <h2 className="font-semibold">{transcript.id}</h2>
            <div className="text-sm text-gray-500">
              {transcript.assembly} &middot; {transcript.type} &middot; {transcript.length} AA
            </div>
          </div>
          <div className="flex gap-2">
            <button onClick={copy}     className="px-3 py-1.5 text-sm border rounded hover:bg-gray-50">Copy</button>
            <button onClick={download} className="px-3 py-1.5 text-sm border rounded hover:bg-gray-50">Download .fasta</button>
            <button onClick={onClose}  className="px-3 py-1.5 text-sm border rounded hover:bg-gray-50">Close</button>
          </div>
        </div>
        <div className="overflow-auto p-4">
          <pre className="font-mono text-xs whitespace-pre-wrap break-all text-gray-700">{fasta}</pre>
        </div>
      </div>
    </div>
  )
}
