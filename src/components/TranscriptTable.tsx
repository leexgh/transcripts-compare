import React, { useState, useCallback, useMemo, useEffect, useRef } from 'react'
import {
  useReactTable, getCoreRowModel, getSortedRowModel, getFilteredRowModel,
  flexRender, type ColumnDef, type SortingState,
} from '@tanstack/react-table'
import { useVirtualizer } from '@tanstack/react-virtual'
import { useNavigate } from 'react-router-dom'
import type { GeneEntry, GeneData, Similarities } from '../lib/types'
import { computeSimilarity } from '../lib/similarity'
import SimilarityBadge from './SimilarityBadge'
import TranscriptDropdown from './TranscriptDropdown'

// Per-row "best match" state (allows dropdown swaps)
interface RowState {
  grch37_enst: string
  grch38_enst: string
  grch37_nm: string
  grch38_nm: string
  similarities: Similarities
}

function buildSims(
  row: GeneEntry,
  bm: { grch37_enst: string; grch38_enst: string; grch37_nm: string; grch38_nm: string },
  seqMap: (id: string, gene: GeneEntry) => string,
): Similarities {
  const s37e = seqMap(bm.grch37_enst, row)
  const s38e = seqMap(bm.grch38_enst, row)
  const s37n = seqMap(bm.grch37_nm, row)
  const s38n = seqMap(bm.grch38_nm, row)
  const mane = seqMap(row.mane_grch38.enst, row)
  return {
    grch37_enst_vs_grch37_refseq: computeSimilarity(s37e, s37n),
    grch38_enst_vs_grch38_refseq: computeSimilarity(s38e, s38n),
    grch37_enst_vs_grch38_enst: computeSimilarity(s37e, s38e),
    grch37_enst_vs_mane: computeSimilarity(s37e, mane),
    grch38_enst_vs_mane: computeSimilarity(s38e, mane),
  }
}

function getSeq(id: string, gene: GeneEntry): string {
  const t = gene.transcripts.find(t => t.id === id)
  return t?.sequence ?? ''
}

interface Props {
  data: GeneData
  genes: GeneEntry[]
  oncokbOnly?: boolean
}

export default function TranscriptTable({ data, genes, oncokbOnly = false }: Props) {
  const navigate = useNavigate()
  const [sorting, setSorting] = useState<SortingState>([])

  // Per-row state (indexed by gene_symbol)
  const [rowStates, setRowStates] = useState<Record<string, RowState>>({})

  // When oncokbOnly toggles, swap row states to oncokb defaults or back to pipeline defaults
  const prevOncokbOnly = useRef(oncokbOnly)
  useEffect(() => {
    if (prevOncokbOnly.current === oncokbOnly) return
    prevOncokbOnly.current = oncokbOnly
    if (oncokbOnly) {
      // Swap to oncokb_best_match for genes that have it
      const newStates: Record<string, RowState> = {}
      for (const gene of genes) {
        if (gene.oncokb_best_match) {
          const bm = gene.oncokb_best_match
          newStates[gene.gene_symbol] = {
            ...bm,
            similarities: buildSims(gene, bm, getSeq),
          }
        }
      }
      setRowStates(newStates)
    } else {
      // Reset to pipeline defaults
      setRowStates({})
    }
  }, [oncokbOnly, genes])

  const getRowState = useCallback((gene: GeneEntry): RowState => {
    return rowStates[gene.gene_symbol] ?? {
      ...gene.best_match,
      similarities: gene.similarities,
    }
  }, [rowStates])

  const updateRow = useCallback((gene: GeneEntry, patch: Partial<RowState>) => {
    setRowStates(prev => {
      const cur = prev[gene.gene_symbol] ?? { ...gene.best_match, similarities: gene.similarities }
      const next = { ...cur, ...patch }
      next.similarities = buildSims(gene, next, (id, g) => getSeq(id, g))
      return { ...prev, [gene.gene_symbol]: next }
    })
  }, [])

  const columns = useMemo<ColumnDef<GeneEntry>[]>(() => [
    {
      accessorKey: 'gene_symbol',
      header: 'Gene',
      size: 90,
      cell: ({ getValue, row }) => (
        <div>
          <span className="font-semibold">{getValue() as string}</span>
          <div className="flex gap-1 mt-0.5">
            {row.original.is_germline && <span className="px-1 py-0.5 bg-purple-100 text-purple-700 text-[10px] rounded">GL</span>}
            {row.original.is_clinical && <span className="px-1 py-0.5 bg-blue-100 text-blue-700 text-[10px] rounded">CL</span>}
          </div>
        </div>
      ),
    },
    {
      id: 'hgnc_entrez',
      header: 'HGNC / Entrez',
      size: 130,
      cell: ({ row }) => {
        const hgnc = row.original.hgnc_id
        const entrez = row.original.entrez_gene_id
        if (!hgnc && !entrez) return <span className="text-gray-300">—</span>
        const hgncNum = hgnc ? hgnc.replace('HGNC:', '') : ''
        const label = [hgncNum ? `HGNC:${hgncNum}` : '', entrez || ''].filter(Boolean).join('/')
        const href = entrez
          ? `https://www.ncbi.nlm.nih.gov/gene/${entrez}`
          : hgnc
            ? `https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/${hgnc}`
            : ''
        return href ? (
          <a href={href} target="_blank" rel="noopener noreferrer"
             className="text-blue-600 hover:underline text-xs">{label}</a>
        ) : <span className="text-xs">{label}</span>
      },
    },
    {
      id: 'grch37_enst',
      header: 'GRCh37 ENST',
      size: 180,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        return (
          <TranscriptDropdown
            current={rs.grch37_enst}
            transcripts={row.original.transcripts}
            type="ensembl" assembly="GRCh37"
            onChange={t => updateRow(row.original, { grch37_enst: t.id })}
          />
        )
      },
    },
    {
      id: 'grch38_enst',
      header: 'GRCh38 ENST',
      size: 180,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        return (
          <TranscriptDropdown
            current={rs.grch38_enst}
            transcripts={row.original.transcripts}
            type="ensembl" assembly="GRCh38"
            onChange={t => updateRow(row.original, { grch38_enst: t.id })}
          />
        )
      },
    },
    {
      id: 'grch37_nm',
      header: 'GRCh37 NM',
      size: 170,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        return (
          <TranscriptDropdown
            current={rs.grch37_nm}
            transcripts={row.original.transcripts}
            type="refseq" assembly="GRCh37"
            onChange={t => updateRow(row.original, { grch37_nm: t.id })}
          />
        )
      },
    },
    {
      id: 'grch38_nm',
      header: 'GRCh38 NM',
      size: 170,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        return (
          <TranscriptDropdown
            current={rs.grch38_nm}
            transcripts={row.original.transcripts}
            type="refseq" assembly="GRCh38"
            onChange={t => updateRow(row.original, { grch38_nm: t.id })}
          />
        )
      },
    },
    {
      id: 'is_mane',
      header: 'is_MANE',
      size: 140,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        const mane = row.original.mane_grch38
        if (!mane.enst) return <span className="text-gray-300">—</span>
        // Only show badge if currently selected GRCh38 ENST matches the MANE ENST
        const maneBase = mane.enst.split('.')[0]
        const selectedBase = rs.grch38_enst.split('.')[0]
        if (maneBase !== selectedBase) return <span className="text-gray-300">—</span>
        const isSelect = mane.status === 'MANE Select'
        return (
          <span className={`px-1.5 py-0.5 rounded text-[11px] ${isSelect ? 'bg-green-100 text-green-800' : 'bg-teal-100 text-teal-800'}`}>
            {isSelect ? 'Select' : 'Plus Clinical'}
          </span>
        )
      },
    },
    // 5 similarity columns
    {
      id: 'sim_37e_37r',
      header: '37E/37R',
      size: 100,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        return (
          <SimilarityBadge
            sim={rs.similarities.grch37_enst_vs_grch37_refseq}
            label="GRCh37 ENST vs GRCh37 RefSeq"
            onClick={() => navigate(`/diff?a=${encodeURIComponent(rs.grch37_enst)}&b=${encodeURIComponent(rs.grch37_nm)}`)}
          />
        )
      },
    },
    {
      id: 'sim_38e_38r',
      header: '38E/38R',
      size: 100,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        return (
          <SimilarityBadge
            sim={rs.similarities.grch38_enst_vs_grch38_refseq}
            label="GRCh38 ENST vs GRCh38 RefSeq"
            onClick={() => navigate(`/diff?a=${encodeURIComponent(rs.grch38_enst)}&b=${encodeURIComponent(rs.grch38_nm)}`)}
          />
        )
      },
    },
    {
      id: 'sim_37e_38e',
      header: '37E/38E',
      size: 100,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        return (
          <SimilarityBadge
            sim={rs.similarities.grch37_enst_vs_grch38_enst}
            label="GRCh37 ENST vs GRCh38 ENST"
            onClick={() => navigate(`/diff?a=${encodeURIComponent(rs.grch37_enst)}&b=${encodeURIComponent(rs.grch38_enst)}`)}
          />
        )
      },
    },
    {
      id: 'sim_37e_mane',
      header: '37E/MANE',
      size: 105,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        return (
          <SimilarityBadge
            sim={rs.similarities.grch37_enst_vs_mane}
            label="GRCh37 ENST vs MANE GRCh38"
            onClick={() => navigate(`/diff?a=${encodeURIComponent(rs.grch37_enst)}&b=${encodeURIComponent(row.original.mane_grch38.enst)}`)}
          />
        )
      },
    },
    {
      id: 'sim_38e_mane',
      header: '38E/MANE',
      size: 105,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        return (
          <SimilarityBadge
            sim={rs.similarities.grch38_enst_vs_mane}
            label="GRCh38 ENST vs MANE GRCh38"
            onClick={() => navigate(`/diff?a=${encodeURIComponent(rs.grch38_enst)}&b=${encodeURIComponent(row.original.mane_grch38.enst)}`)}
          />
        )
      },
    },
    {
      accessorKey: 'curated_note',
      header: 'Note',
      size: 160,
      cell: ({ getValue }) => <span className="text-xs text-gray-500">{getValue() as string}</span>,
    },
  ], [getRowState, updateRow, navigate])

  const table = useReactTable({
    data: genes,
    columns,
    state: { sorting },
    onSortingChange: setSorting,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
  })

  const { rows } = table.getRowModel()
  const parentRef = React.useRef<HTMLDivElement>(null)

  const virtualizer = useVirtualizer({
    count: rows.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 48,
    overscan: 10,
  })

  const totalSize = virtualizer.getTotalSize()
  const virtualRows = virtualizer.getVirtualItems()

  return (
    <div className="flex flex-col flex-1 overflow-hidden">
      <div ref={parentRef} className="overflow-auto flex-1">
        <table className="min-w-full text-sm border-collapse">
          <thead className="bg-gray-50 sticky top-0 z-10">
            {table.getHeaderGroups().map(hg => (
              <tr key={hg.id}>
                {hg.headers.map(h => (
                  <th
                    key={h.id}
                    style={{ width: h.getSize() }}
                    className="text-left px-3 py-2 text-xs font-medium text-gray-600 border-b whitespace-nowrap cursor-pointer select-none hover:bg-gray-100"
                    onClick={h.column.getToggleSortingHandler()}
                  >
                    {flexRender(h.column.columnDef.header, h.getContext())}
                    {{ asc: ' ↑', desc: ' ↓' }[h.column.getIsSorted() as string] ?? ''}
                  </th>
                ))}
              </tr>
            ))}
          </thead>
          <tbody>
            <tr><td style={{ height: virtualRows[0]?.start ?? 0 }} /></tr>
            {virtualRows.map(vr => {
              const row = rows[vr.index]
              return (
                <tr key={row.id} data-index={vr.index} ref={virtualizer.measureElement} className="border-b hover:bg-gray-50">
                  {row.getVisibleCells().map(cell => (
                    <td key={cell.id} className="px-3 py-2 align-top" style={{ width: cell.column.getSize() }}>
                      {flexRender(cell.column.columnDef.cell, cell.getContext())}
                    </td>
                  ))}
                </tr>
              )
            })}
            <tr><td style={{ height: totalSize - (virtualRows[virtualRows.length - 1]?.end ?? totalSize) }} /></tr>
          </tbody>
        </table>
      </div>

    </div>
  )
}
