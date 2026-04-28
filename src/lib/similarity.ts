import type { SimilarityResult } from './types'

/**
 * Compute sequence similarity using a longest-common-subsequence approach.
 * Uses the same formula as Python's difflib.SequenceMatcher.ratio():
 *   2 * lcs / (len1 + len2)
 *
 * The old same-length char-by-char shortcut has been removed. Same-length
 * isoforms can have compensating insertions+deletions; char-by-char treats
 * the entire shifted region as mismatches, producing wildly wrong scores.
 * LCS handles all cases correctly and consistently with precomputed values.
 */
export function computeSimilarity(seq1: string, seq2: string): SimilarityResult {
  if (!seq1 || !seq2) return { pct: null, diff_count: 0 }
  if (seq1 === seq2)  return { pct: 100,  diff_count: 0 }

  const lcsLen     = lcs(seq1, seq2)
  const pct        = Math.round((2 * lcsLen / (seq1.length + seq2.length)) * 10000) / 100
  const diff_count = seq1.length + seq2.length - 2 * lcsLen
  return { pct, diff_count }
}

/** Length of longest common subsequence via space-optimised DP. */
function lcs(a: string, b: string): number {
  const m = a.length, n = b.length
  let prev = new Uint16Array(n + 1)
  let curr = new Uint16Array(n + 1)
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      curr[j] = a[i - 1] === b[j - 1] ? prev[j - 1] + 1 : Math.max(prev[j], curr[j - 1])
    }
    ;[prev, curr] = [curr, prev]
    curr.fill(0)
  }
  return prev[n]
}
