import { useState, useCallback } from 'react'
import type { Transcript } from '../lib/types'

export function useCompareSelection() {
  const [selected, setSelected] = useState<Transcript[]>([])

  const toggle = useCallback((t: Transcript) => {
    setSelected(prev => {
      const idx = prev.findIndex(x => x.id === t.id)
      if (idx >= 0) return prev.filter((_, i) => i !== idx)
      return [...prev, t]
    })
  }, [])

  const clear = useCallback(() => setSelected([]), [])

  const isSelected = useCallback((id: string) =>
    selected.some(t => t.id === id), [selected])

  return { selected, toggle, clear, isSelected }
}
