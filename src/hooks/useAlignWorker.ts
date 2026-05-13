import { useRef, useState, useCallback, useEffect } from 'react'
import type { AlignmentResult } from '../lib/alignment'
import { LONG_SEQ_THRESHOLD } from '../lib/alignment'

export type AlignStatus = 'idle' | 'running' | 'done' | 'error'

export function useAlignWorker(seq1: string, seq2: string) {
  const isLong = seq1.length > LONG_SEQ_THRESHOLD || seq2.length > LONG_SEQ_THRESHOLD
  const [status, setStatus] = useState<AlignStatus>('idle')
  const [result, setResult] = useState<AlignmentResult | null>(null)
  const workerRef = useRef<Worker | null>(null)
  const seqsRef = useRef({ seq1, seq2 })
  seqsRef.current = { seq1, seq2 }

  const run = useCallback(() => {
    workerRef.current?.terminate()
    setStatus('running')
    setResult(null)
    const { seq1: s1, seq2: s2 } = seqsRef.current
    const worker = new Worker(
      new URL('../workers/alignWorker.ts', import.meta.url),
      { type: 'module' },
    )
    workerRef.current = worker
    worker.onmessage = (e: MessageEvent<AlignmentResult>) => {
      setResult(e.data)
      setStatus('done')
    }
    worker.onerror = () => setStatus('error')
    worker.postMessage({ seq1: s1, seq2: s2 })
  }, [])

  useEffect(() => {
    if (seq1 && seq2 && !isLong) {
      run()
    } else {
      setStatus('idle')
      setResult(null)
    }
    return () => workerRef.current?.terminate()
  }, [seq1, seq2, isLong, run])

  return { result, status, isLong, run }
}
