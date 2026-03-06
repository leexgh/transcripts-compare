/**
 * Needleman-Wunsch global alignment for protein sequences.
 * Falls back to char-by-char for same-length sequences.
 */

export interface AlignmentResult {
  seq1Aligned: string
  seq2Aligned: string
  score:       number
  pct:         number
  diffCount:   number
}

const GAP_OPEN = -10
const GAP_EXT  = -1

// Simplified scoring: match=1, mismatch=-1, gap=GAP_OPEN
function score(a: string, b: string): number {
  if (a === '-' || b === '-') return GAP_OPEN
  return a === b ? 1 : -1
}

export function align(seq1: string, seq2: string): AlignmentResult {
  if (seq1.length === seq2.length) {
    // Simple char-by-char
    let diffs = 0
    const parts1: string[] = []
    const parts2: string[] = []
    for (let i = 0; i < seq1.length; i++) {
      parts1.push(seq1[i])
      parts2.push(seq2[i])
      if (seq1[i] !== seq2[i]) diffs++
    }
    const pct = seq1.length > 0
      ? Math.round(((seq1.length - diffs) / seq1.length) * 10000) / 100
      : 100
    return { seq1Aligned: seq1, seq2Aligned: seq2, score: seq1.length - diffs, pct, diffCount: diffs }
  }

  // Needleman-Wunsch for different lengths
  const n = seq1.length
  const m = seq2.length

  // Cap at 2000 AA for performance
  if (n > 2000 || m > 2000) {
    return alignTruncated(seq1.slice(0, 2000), seq2.slice(0, 2000))
  }

  const dp: number[][] = Array.from({ length: n + 1 }, () => new Array(m + 1).fill(0))
  // Init
  for (let i = 0; i <= n; i++) dp[i][0] = i * GAP_OPEN
  for (let j = 0; j <= m; j++) dp[0][j] = j * GAP_OPEN

  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      const match = dp[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1])
      const del   = dp[i - 1][j]     + GAP_OPEN
      const ins   = dp[i][j - 1]     + GAP_OPEN
      dp[i][j] = Math.max(match, del, ins)
    }
  }

  // Traceback
  let i = n, j = m
  let a1 = '', a2 = ''
  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 && dp[i][j] === dp[i-1][j-1] + score(seq1[i-1], seq2[j-1])) {
      a1 = seq1[i-1] + a1
      a2 = seq2[j-1] + a2
      i--; j--
    } else if (i > 0 && dp[i][j] === dp[i-1][j] + GAP_OPEN) {
      a1 = seq1[i-1] + a1
      a2 = '-' + a2
      i--
    } else {
      a1 = '-' + a1
      a2 = seq2[j-1] + a2
      j--
    }
  }

  let diffs = 0
  let totalCols = a1.length
  for (let k = 0; k < totalCols; k++) {
    if (a1[k] !== a2[k]) diffs++
  }
  const pct = totalCols > 0
    ? Math.round(((totalCols - diffs) / totalCols) * 10000) / 100
    : 100

  return { seq1Aligned: a1, seq2Aligned: a2, score: dp[n][m], pct, diffCount: diffs }
}

function alignTruncated(s1: string, s2: string): AlignmentResult {
  // Simple alignment for truncated sequences
  const result = align(s1, s2)
  return result
}

export interface DiffBlock {
  type:     'match' | 'sub' | 'ins' | 'del'
  pos1:     number
  pos2:     number
  seq1chars: string
  seq2chars: string
}

export function buildDiffBlocks(seq1Aligned: string, seq2Aligned: string): DiffBlock[] {
  const blocks: DiffBlock[] = []
  let pos1 = 0, pos2 = 0
  const n = seq1Aligned.length

  for (let i = 0; i < n; i++) {
    const c1 = seq1Aligned[i]
    const c2 = seq2Aligned[i]
    let type: DiffBlock['type']
    if (c1 === c2)             type = 'match'
    else if (c1 === '-')       type = 'ins'
    else if (c2 === '-')       type = 'del'
    else                       type = 'sub'

    // Merge consecutive same-type blocks
    const last = blocks[blocks.length - 1]
    if (last && last.type === type && type !== 'match') {
      last.seq1chars += c1
      last.seq2chars += c2
    } else if (type !== 'match') {
      blocks.push({ type, pos1, pos2, seq1chars: c1, seq2chars: c2 })
    }

    if (c1 !== '-') pos1++
    if (c2 !== '-') pos2++
  }

  return blocks
}
