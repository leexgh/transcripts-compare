import type { SimilarityResult } from './types'

/**
 * Compute sequence similarity using a longest-common-subsequence approach.
 * Correctly handles insertions and deletions (not just char-by-char).
 * Uses the same formula as Python's difflib.SequenceMatcher.ratio():
 *   2 * matches / (len1 + len2)
 */
export function computeSimilarity(seq1: string, seq2: string): SimilarityResult {
  if (!seq1 || !seq2) return { pct: null, diff_count: 0 }

  if (seq1 === seq2) return { pct: 100, diff_count: 0 }

  // For same-length sequences: fast char-by-char
  if (seq1.length === seq2.length) {
    let mismatches = 0
    for (let i = 0; i < seq1.length; i++) {
      if (seq1[i] !== seq2[i]) mismatches++
    }
    const pct = Math.round(((seq1.length - mismatches) / seq1.length) * 10000) / 100
    return { pct, diff_count: mismatches }
  }

  // For different-length sequences: LCS-based (mirrors difflib ratio)
  const lcsLen = lcs(seq1, seq2)
  const pct = Math.round((2 * lcsLen / (seq1.length + seq2.length)) * 10000) / 100
  // diff_count = characters not in LCS (indels + subs)
  const diff_count = seq1.length + seq2.length - 2 * lcsLen
  return { pct, diff_count }
}

/** Length of longest common subsequence (capped at 1000 chars for performance) */
function lcs(a: string, b: string): number {
  const maxLen = 1000
  const s1 = a.length > maxLen ? a.slice(0, maxLen) : a
  const s2 = b.length > maxLen ? b.slice(0, maxLen) : b
  const m = s1.length, n = s2.length
  // Space-optimised DP (two rows)
  let prev = new Uint16Array(n + 1)
  let curr = new Uint16Array(n + 1)
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      curr[j] = s1[i - 1] === s2[j - 1] ? prev[j - 1] + 1 : Math.max(prev[j], curr[j - 1])
    }
    ;[prev, curr] = [curr, prev]
    curr.fill(0)
  }
  // Scale back if we truncated
  const scale = Math.min(a.length, maxLen) === a.length ? 1 : a.length / maxLen
  return Math.round(prev[n] * scale)
}
