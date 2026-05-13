import { align } from '../lib/alignment'
import type { AlignmentResult } from '../lib/alignment'

self.onmessage = (e: MessageEvent<{ seq1: string; seq2: string }>) => {
  const result: AlignmentResult = align(e.data.seq1, e.data.seq2)
  self.postMessage(result)
}
