import React from 'react'
import type { FilterState, AssemblyFilter } from '../lib/types'

interface Props {
  filters:   FilterState
  onChange:  (f: FilterState) => void
  geneCount: number
  total:     number
}

export default function FilterBar({ filters, onChange, geneCount, total }: Props) {
  const set = (patch: Partial<FilterState>) => onChange({ ...filters, ...patch })

  return (
    <div className="bg-white border-b px-4 py-3 flex flex-wrap gap-3 items-center">
      {/* Search */}
      <input
        type="search"
        placeholder="Search gene / transcript ID…"
        value={filters.search}
        onChange={e => set({ search: e.target.value })}
        className="border rounded px-3 py-1.5 text-sm w-56 focus:outline-none focus:ring-2 focus:ring-blue-300"
      />

      {/* Assembly toggle */}
      <div className="flex border rounded overflow-hidden text-sm">
        {(['both', 'GRCh37', 'GRCh38'] as AssemblyFilter[]).map(a => (
          <button
            key={a}
            onClick={() => set({ assembly: a })}
            className={`px-3 py-1.5 ${filters.assembly === a ? 'bg-blue-600 text-white' : 'hover:bg-gray-50'}`}
          >
            {a === 'both' ? 'Both' : a}
          </button>
        ))}
      </div>

      {/* Checkboxes */}
      {[
        { key: 'clinicalOnly',  label: 'Clinical only'  },
        { key: 'germlineOnly',  label: 'Germline only'  },
        { key: 'maneOnly',      label: 'MANE only'      },
        { key: 'mismatchOnly',  label: 'Mismatches only' },
      ].map(({ key, label }) => (
        <label key={key} className="flex items-center gap-1.5 text-sm cursor-pointer">
          <input
            type="checkbox"
            checked={filters[key as keyof FilterState] as boolean}
            onChange={e => set({ [key]: e.target.checked })}
            className="rounded"
          />
          {label}
        </label>
      ))}

      <span className="ml-auto text-sm text-gray-500">
        {geneCount} / {total} genes
      </span>
    </div>
  )
}
