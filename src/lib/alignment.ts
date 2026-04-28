/**
 * Needleman-Wunsch global alignment for protein sequences.
 *
 * Always uses NW — the same-length char-by-char shortcut was removed because
 * isoforms can have the same total length while containing a compensating
 * insertion+deletion. In that case char-by-char mislabels the entire shifted
 * region as substitutions instead of one ins + one del, producing a wildly
 * wrong identity score.
 *
 * Percent identity uses the same formula as Python difflib / computeSimilarity:
 *   2 * matches / (len1 + len2)
 * so the score shown in the diff viewer is consistent with the table.
 */

export interface AlignmentResult {
  seq1Aligned: string
  seq2Aligned: string
  score:       number
  pct:         number
  diffCount:   number
}

const GAP_OPEN = -10

export function align(seq1: string, seq2: string): AlignmentResult {
  const n = seq1.length
  const m = seq2.length

  if (n === 0 || m === 0) {
    return { seq1Aligned: seq1, seq2Aligned: seq2, score: 0, pct: 0, diffCount: n + m }
  }

  // Cap at 2000 AA each for performance (O(n·m) memory + time).
  if (n > 2000 || m > 2000) {
    return align(seq1.slice(0, 2000), seq2.slice(0, 2000))
  }

  // ── Needleman-Wunsch DP (linear gap penalty) ──────────────────────────────
  // Use a flat Int32Array to avoid the overhead of a jagged array.
  const W = m + 1
  const dp = new Int32Array((n + 1) * W)

  for (let i = 0; i <= n; i++) dp[i * W] = i * GAP_OPEN
  for (let j = 0; j <= m; j++) dp[j]     = j * GAP_OPEN

  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      const sub   = dp[(i - 1) * W + (j - 1)] + (seq1[i - 1] === seq2[j - 1] ? 1 : -1)
      const del_  = dp[(i - 1) * W + j]       + GAP_OPEN
      const ins   = dp[i       * W + (j - 1)] + GAP_OPEN
      dp[i * W + j] = sub > del_ ? (sub > ins ? sub : ins) : (del_ > ins ? del_ : ins)
    }
  }

  // ── Traceback ─────────────────────────────────────────────────────────────
  let i = n, j = m
  let a1 = '', a2 = ''
  while (i > 0 || j > 0) {
    if (
      i > 0 && j > 0 &&
      dp[i * W + j] === dp[(i - 1) * W + (j - 1)] + (seq1[i - 1] === seq2[j - 1] ? 1 : -1)
    ) {
      a1 = seq1[i - 1] + a1
      a2 = seq2[j - 1] + a2
      i--; j--
    } else if (i > 0 && dp[i * W + j] === dp[(i - 1) * W + j] + GAP_OPEN) {
      a1 = seq1[i - 1] + a1
      a2 = '-' + a2
      i--
    } else {
      a1 = '-' + a1
      a2 = seq2[j - 1] + a2
      j--
    }
  }

  // ── Metrics (consistent with Python difflib / computeSimilarity) ──────────
  // Count only non-gap identical columns; divide by total original residues.
  let matches = 0
  for (let k = 0; k < a1.length; k++) {
    if (a1[k] !== '-' && a1[k] === a2[k]) matches++
  }
  const pct       = Math.round((2 * matches / (n + m)) * 10000) / 100
  const diffCount = n + m - 2 * matches

  return { seq1Aligned: a1, seq2Aligned: a2, score: dp[n * W + m], pct, diffCount }
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
