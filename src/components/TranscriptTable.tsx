import React, { useState, useCallback, useMemo, useRef } from 'react'
import {
  useReactTable, getCoreRowModel, getSortedRowModel, getFilteredRowModel,
  flexRender, type ColumnDef, type SortingState, type SortingFn,
} from '@tanstack/react-table'
import { useVirtualizer } from '@tanstack/react-virtual'
import { useSearchParams } from 'react-router-dom'
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
  const s37e    = seqMap(bm.grch37_enst, row)
  const s38e    = seqMap(bm.grch38_enst, row)
  const s37n    = seqMap(bm.grch37_nm, row)
  const s38n    = seqMap(bm.grch38_nm, row)
  const sMane   = seqMap(row.mane_grch38.enst, row)
  const sManeNm = seqMap(row.mane_grch38.nm, row)
  return {
    grch37_enst_vs_grch37_refseq:   computeSimilarity(s37e, s37n),
    grch38_enst_vs_grch38_refseq:   computeSimilarity(s38e, s38n),
    grch37_enst_vs_grch38_enst:     computeSimilarity(s37e, s38e),
    grch37_enst_vs_mane:            computeSimilarity(s37e, sMane),
    grch38_enst_vs_mane:            computeSimilarity(s38e, sMane),
    grch37_nm_vs_mane:              computeSimilarity(s37n, sMane),
    grch37_enst_vs_grch38_refseq:   computeSimilarity(s37e, s38n),
    grch37_enst_vs_mane_nm:         computeSimilarity(s37e, sManeNm),
    grch38_enst_vs_grch37_refseq:   computeSimilarity(s38e, s37n),
    grch38_enst_vs_mane_nm:         computeSimilarity(s38e, sManeNm),
    grch37_refseq_vs_grch38_refseq: computeSimilarity(s37n, s38n),
    grch37_refseq_vs_mane_nm:       computeSimilarity(s37n, sManeNm),
    grch38_refseq_vs_mane:          computeSimilarity(s38n, sMane),
    grch38_refseq_vs_mane_nm:       computeSimilarity(s38n, sManeNm),
    mane_vs_mane_nm:                computeSimilarity(sMane, sManeNm),
  }
}

function getSeq(id: string, gene: GeneEntry): string {
  const t = gene.transcripts.find(t => t.id === id)
  return t?.sequence ?? ''
}

// Custom sorting function: NA (-1) always sorts last regardless of direction.
// TanStack Table inverts the comparator result when sorting descending, so we
// read the current sort direction (via ref) and return a pre-inverted value for NA rows.
function makeSimSortingFn(sortingRef: React.RefObject<SortingState>): SortingFn<GeneEntry> {
  return (rowA, rowB, columnId) => {
    const a = rowA.getValue<number>(columnId)
    const b = rowB.getValue<number>(columnId)
    const aIsNA = a === -1
    const bIsNA = b === -1
    if (aIsNA && bIsNA) return 0
    if (aIsNA || bIsNA) {
      const desc = (sortingRef.current ?? []).find(s => s.id === columnId)?.desc ?? false
      if (aIsNA) return desc ? -1 : 1   // post-inversion: always positive → a last
      return desc ? 1 : -1               // post-inversion: always negative → b last
    }
    return a - b
  }
}

type CompareOption = '37ENST' | '38ENST' | '37RefSeq' | '38RefSeq' | '38MANE(ENST)' | '38MANE(RefSeq)'

const COMPARE_OPTIONS: CompareOption[] = [
  '37ENST', '38ENST', '37RefSeq', '38RefSeq', '38MANE(ENST)', '38MANE(RefSeq)',
]

// Short abbreviations used in URL param to keep it clean
const OPTION_ABBREV: Record<CompareOption, string> = {
  '37ENST': '37e', '38ENST': '38e',
  '37RefSeq': '37r', '38RefSeq': '38r',
  '38MANE(ENST)': '38me', '38MANE(RefSeq)': '38mr',
}
const ABBREV_OPTION: Record<string, CompareOption> = Object.fromEntries(
  Object.entries(OPTION_ABBREV).map(([k, v]) => [v, k as CompareOption])
)

const DEFAULT_COLS_PARAM = '37e~37r,38e~38r,37e~38e,37e~38me,38e~38me,37r~38me'

interface ActiveCol {
  id: string          // stable TanStack column id, e.g. 'sim_37e~37r'
  a: CompareOption
  b: CompareOption
}

function parseCols(param: string | null): ActiveCol[] {
  const raw = param || DEFAULT_COLS_PARAM
  const counts: Record<string, number> = {}
  return raw.split(',').flatMap(part => {
    const [abbrevA, abbrevB] = part.split('~')
    const a = ABBREV_OPTION[abbrevA]
    const b = ABBREV_OPTION[abbrevB]
    if (!a || !b) return []
    counts[part] = (counts[part] ?? 0) + 1
    const id = counts[part] === 1 ? `sim_${part}` : `sim_${part}_${counts[part]}`
    return [{ id, a, b }]
  })
}

// Map every (A, B) dropdown pair to a precomputed Similarities key (both directions).
const PRECOMPUTED_SIM_KEY: Record<string, keyof Similarities> = {
  // original 6 pairs
  '37ENST|37RefSeq':         'grch37_enst_vs_grch37_refseq',
  '37RefSeq|37ENST':         'grch37_enst_vs_grch37_refseq',
  '38ENST|38RefSeq':         'grch38_enst_vs_grch38_refseq',
  '38RefSeq|38ENST':         'grch38_enst_vs_grch38_refseq',
  '37ENST|38ENST':           'grch37_enst_vs_grch38_enst',
  '38ENST|37ENST':           'grch37_enst_vs_grch38_enst',
  '37ENST|38MANE(ENST)':     'grch37_enst_vs_mane',
  '38MANE(ENST)|37ENST':     'grch37_enst_vs_mane',
  '38ENST|38MANE(ENST)':     'grch38_enst_vs_mane',
  '38MANE(ENST)|38ENST':     'grch38_enst_vs_mane',
  '37RefSeq|38MANE(ENST)':   'grch37_nm_vs_mane',
  '38MANE(ENST)|37RefSeq':   'grch37_nm_vs_mane',
  // new 9 pairs
  '37ENST|38RefSeq':         'grch37_enst_vs_grch38_refseq',
  '38RefSeq|37ENST':         'grch37_enst_vs_grch38_refseq',
  '37ENST|38MANE(RefSeq)':   'grch37_enst_vs_mane_nm',
  '38MANE(RefSeq)|37ENST':   'grch37_enst_vs_mane_nm',
  '38ENST|37RefSeq':         'grch38_enst_vs_grch37_refseq',
  '37RefSeq|38ENST':         'grch38_enst_vs_grch37_refseq',
  '38ENST|38MANE(RefSeq)':   'grch38_enst_vs_mane_nm',
  '38MANE(RefSeq)|38ENST':   'grch38_enst_vs_mane_nm',
  '37RefSeq|38RefSeq':       'grch37_refseq_vs_grch38_refseq',
  '38RefSeq|37RefSeq':       'grch37_refseq_vs_grch38_refseq',
  '37RefSeq|38MANE(RefSeq)': 'grch37_refseq_vs_mane_nm',
  '38MANE(RefSeq)|37RefSeq': 'grch37_refseq_vs_mane_nm',
  '38RefSeq|38MANE(ENST)':   'grch38_refseq_vs_mane',
  '38MANE(ENST)|38RefSeq':   'grch38_refseq_vs_mane',
  '38RefSeq|38MANE(RefSeq)': 'grch38_refseq_vs_mane_nm',
  '38MANE(RefSeq)|38RefSeq': 'grch38_refseq_vs_mane_nm',
  '38MANE(ENST)|38MANE(RefSeq)': 'mane_vs_mane_nm',
  '38MANE(RefSeq)|38MANE(ENST)': 'mane_vs_mane_nm',
}

function getCompareId(option: CompareOption, rs: RowState, gene: GeneEntry): string {
  switch (option) {
    case '37ENST':         return rs.grch37_enst
    case '38ENST':         return rs.grch38_enst
    case '37RefSeq':       return rs.grch37_nm
    case '38RefSeq':       return rs.grch38_nm
    case '38MANE(ENST)':   return gene.mane_grch38.enst
    case '38MANE(RefSeq)': return gene.mane_grch38.nm
  }
}

interface Props {
  data: GeneData
  genes: GeneEntry[]
  sorting: SortingState
  onSortingChange: (updater: SortingState | ((prev: SortingState) => SortingState)) => void
}

function openDiff(a: string, b: string) {
  window.open(
    `${window.location.href.split('#')[0]}#/diff?a=${encodeURIComponent(a)}&b=${encodeURIComponent(b)}`,
    '_blank',
  )
}

export default function TranscriptTable({ data, genes, sorting, onSortingChange }: Props) {
  const [searchParams, setSearchParams] = useSearchParams()

  const activeCols = useMemo(() => parseCols(searchParams.get('cols')), [searchParams])

  // Update the cols URL param, optionally patching sort in the same call
  const setColsParam = (
    newCols: Array<{ a: CompareOption; b: CompareOption }>,
    patchSort?: (p: URLSearchParams) => void,
  ) => {
    setSearchParams(prev => {
      const p = new URLSearchParams(prev)
      const serialized = newCols.map(c => `${OPTION_ABBREV[c.a]}~${OPTION_ABBREV[c.b]}`).join(',')
      if (serialized === DEFAULT_COLS_PARAM) p.delete('cols')
      else p.set('cols', serialized)
      patchSort?.(p)
      return p
    }, { replace: true })
  }

  const removeCol = (colIndex: number) => {
    const removedId = activeCols[colIndex].id
    setColsParam(
      activeCols.filter((_, i) => i !== colIndex),
      p => {
        const parts = (p.get('sort') ?? '').split(',').filter(s => s && !s.startsWith(removedId + ':'))
        if (parts.length === 0) p.delete('sort'); else p.set('sort', parts.join(','))
      },
    )
  }

  const updateColAB = (colIndex: number, patch: { a?: CompareOption; b?: CompareOption }) => {
    const oldCol = activeCols[colIndex]
    const oldId = oldCol.id
    const updated = activeCols.map((c, i) => i === colIndex ? { ...c, ...patch } : c)
    // Recompute the new id the same way parseCols would
    const newId = parseCols(updated.map(c => `${OPTION_ABBREV[c.a]}~${OPTION_ABBREV[c.b]}`).join(','))[colIndex].id
    setColsParam(
      updated,
      p => {
        if (oldId !== newId) {
          const parts = (p.get('sort') ?? '').split(',').filter(Boolean)
            .map(s => s.startsWith(oldId + ':') ? s.replace(oldId + ':', newId + ':') : s)
          if (parts.length === 0) p.delete('sort'); else p.set('sort', parts.join(','))
        }
      },
    )
  }

  const addCol = () => setColsParam([...activeCols, { a: '37ENST', b: '38ENST' }])

  const sortingRef = useRef<SortingState>(sorting)
  sortingRef.current = sorting
  const simSortingFn = useMemo(() => makeSimSortingFn(sortingRef), [])

  // Per-row state (indexed by gene_symbol)
  const [rowStates, setRowStates] = useState<Record<string, RowState>>({})

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

  // Precompute similarities for all compare columns so sorting is O(1) per comparison.
  // Uses precomputed gene.similarities when the pair matches a known key to avoid expensive LCS.
  const simMaps = useMemo(() => {
    const maps: Record<string, Map<string, number>> = {}
    for (const col of activeCols) {
      const preKey = PRECOMPUTED_SIM_KEY[`${col.a}|${col.b}`]
      const map = new Map<string, number>()
      for (const gene of genes) {
        const rs = getRowState(gene)
        const preSim = preKey ? rs.similarities[preKey] : undefined
        if (preSim !== undefined) {
          map.set(gene.gene_symbol, preSim.pct ?? -1)
        } else {
          const idA = getCompareId(col.a, rs, gene)
          const idB = getCompareId(col.b, rs, gene)
          if (!idA || !idB) { map.set(gene.gene_symbol, -1); continue }
          map.set(gene.gene_symbol, computeSimilarity(getSeq(idA, gene), getSeq(idB, gene)).pct ?? -1)
        }
      }
      maps[col.id] = map
    }
    return maps
  }, [genes, activeCols, getRowState])

  const columns = useMemo<ColumnDef<GeneEntry>[]>(() => [
    {
      accessorKey: 'gene_symbol',
      header: 'Gene',
      size: 90,
      cell: ({ getValue, row }) => (
        <div>
          <span className="font-semibold">{getValue() as string}</span>
          <div className="flex gap-1 mt-0.5">
            {row.original.is_germline && <span className="px-1 py-0.5 bg-purple-100 text-purple-700 text-[10px] rounded">Germline</span>}
            {row.original.is_clinical && <span className="px-1 py-0.5 bg-blue-100 text-blue-700 text-[10px] rounded">Clinical</span>}
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
      header: '37ENST',
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
      header: '38ENST',
      size: 180,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        const mane = row.original.mane_grch38
        const maneEnstMatch = mane.enst && mane.enst.split('.')[0] === rs.grch38_enst.split('.')[0]
        return (
          <div>
            <TranscriptDropdown
              current={rs.grch38_enst}
              transcripts={row.original.transcripts}
              type="ensembl" assembly="GRCh38"
              onChange={t => updateRow(row.original, { grch38_enst: t.id })}
            />
            {maneEnstMatch && (
              <span className="px-1 py-0.5 rounded text-[10px] bg-pink-100 text-pink-700">
                MANE
              </span>
            )}
          </div>
        )
      },
    },
    {
      id: 'grch37_nm',
      header: '37RefSeq',
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
      header: '38RefSeq',
      size: 170,
      cell: ({ row }) => {
        const rs = getRowState(row.original)
        const mane = row.original.mane_grch38
        const maneNmMatch = mane.nm && mane.nm.split('.')[0] === rs.grch38_nm.split('.')[0]
        return (
          <div>
            <TranscriptDropdown
              current={rs.grch38_nm}
              transcripts={row.original.transcripts}
              type="refseq" assembly="GRCh38"
              onChange={t => updateRow(row.original, { grch38_nm: t.id })}
            />
            {maneNmMatch && (
              <span className="px-1 py-0.5 rounded text-[10px] bg-pink-100 text-pink-700">MANE</span>
            )}
          </div>
        )
      },
    },
    // similarity columns — each column has configurable dropdowns and a close button
    ...activeCols.map((col, colIndex): ColumnDef<GeneEntry> => {
      const simMap = simMaps[col.id]
      return {
        id: col.id,
        size: 165,
        header: () => (
          <div className="flex flex-col gap-0.5 font-normal">
            <div className="flex items-center justify-between gap-0.5">
              <select
                value={col.a}
                onClick={e => e.stopPropagation()}
                onChange={e => updateColAB(colIndex, { a: e.target.value as CompareOption })}
                className="text-[10px] border border-gray-300 rounded px-1 py-0.5 bg-white cursor-pointer flex-1 min-w-0"
              >
                {COMPARE_OPTIONS.map(o => <option key={o} value={o}>{o}</option>)}
              </select>
              <button
                onClick={e => { e.stopPropagation(); removeCol(colIndex) }}
                title="Remove column"
                className="text-gray-300 hover:text-red-500 leading-none text-base px-0.5 shrink-0"
              >×</button>
            </div>
            <span className="text-[10px] text-gray-400 text-center leading-none">vs</span>
            <select
              value={col.b}
              onClick={e => e.stopPropagation()}
              onChange={e => updateColAB(colIndex, { b: e.target.value as CompareOption })}
              className="text-[10px] border border-gray-300 rounded px-1 py-0.5 bg-white cursor-pointer w-full"
            >
              {COMPARE_OPTIONS.map(o => <option key={o} value={o}>{o}</option>)}
            </select>
          </div>
        ),
        accessorFn: (row: GeneEntry) => simMap?.get(row.gene_symbol) ?? -1,
        sortingFn: simSortingFn,
        cell: ({ row }) => {
          const rs = getRowState(row.original)
          const preKey = PRECOMPUTED_SIM_KEY[`${col.a}|${col.b}`]
          const preSim = preKey ? rs.similarities[preKey] : undefined
          const idA = getCompareId(col.a, rs, row.original)
          const idB = getCompareId(col.b, rs, row.original)
          const sim = preSim !== undefined
            ? preSim
            : (idA && idB ? computeSimilarity(getSeq(idA, row.original), getSeq(idB, row.original)) : null)
          if (!sim) return <span className="text-gray-300">—</span>
          return (
            <SimilarityBadge
              sim={sim}
              label={`${col.a} vs ${col.b}`}
              onClick={() => openDiff(idA, idB)}
            />
          )
        },
      }
    }),
    {
      accessorKey: 'curated_note',
      header: 'Note',
      size: 160,
      cell: ({ getValue }) => <span className="text-xs text-gray-500">{getValue() as string}</span>,
    },
    // "Add column" button as a virtual last column
    {
      id: '__add_col',
      size: 44,
      enableSorting: false,
      header: () => (
        <button
          onClick={e => { e.stopPropagation(); addCol() }}
          title="Add comparison column"
          className="w-7 h-7 rounded-full border-2 border-dashed border-gray-300 hover:border-blue-400 text-gray-400 hover:text-blue-500 flex items-center justify-center text-lg leading-none"
        >+</button>
      ),
      cell: () => null,
    },
  ], [getRowState, updateRow, activeCols, simMaps, simSortingFn, removeCol, updateColAB, addCol])

  const table = useReactTable({
    data: genes,
    columns,
    state: { sorting },
    onSortingChange,
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

  function exportTsv() {
    const esc = (s: string | null | undefined) => String(s ?? '').replace(/\t/g, ' ').replace(/\n/g, ' ')

    const fmtTx = (id: string, gene: GeneEntry): string => {
      if (!id) return ''
      const t = gene.transcripts.find(tx => tx.id === id)
      if (t?.source) return `${id}(${t.source})`
      return id
    }

    const headers: string[] = [
      'Gene', 'Is_Clinical', 'Is_Germline', 'HGNC_ID', 'Entrez_ID',
      '37ENST', '38ENST', '37RefSeq', '38RefSeq', '38MANE_ENST', '38MANE_RefSeq',
    ]
    for (const col of activeCols) {
      headers.push(`${col.a} vs ${col.b} Score`, `${col.a} vs ${col.b} Viz_Link`)
    }
    headers.push('Note')

    const baseUrl = window.location.href.split('#')[0]

    const tsvRows = rows.map(row => {
      const gene = row.original
      const rs = getRowState(gene)

      const cells: string[] = [
        esc(gene.gene_symbol),
        gene.is_clinical ? 'Yes' : 'No',
        gene.is_germline ? 'Yes' : 'No',
        esc(gene.hgnc_id),
        esc(gene.entrez_gene_id),
        esc(fmtTx(rs.grch37_enst, gene)),
        esc(fmtTx(rs.grch38_enst, gene)),
        esc(fmtTx(rs.grch37_nm, gene)),
        esc(fmtTx(rs.grch38_nm, gene)),
        esc(fmtTx(gene.mane_grch38.enst, gene)),
        esc(fmtTx(gene.mane_grch38.nm, gene)),
      ]

      for (const col of activeCols) {
        const preKey = PRECOMPUTED_SIM_KEY[`${col.a}|${col.b}`]
        const preSim = preKey ? rs.similarities[preKey] : undefined
        const idA = getCompareId(col.a, rs, gene)
        const idB = getCompareId(col.b, rs, gene)
        const sim = preSim !== undefined
          ? preSim
          : (idA && idB ? computeSimilarity(getSeq(idA, gene), getSeq(idB, gene)) : null)

        const score = !sim ? '' : sim.pct === null ? 'N/A' : `${sim.pct}%`
        const link = (idA && idB)
          ? `${baseUrl}#/diff?a=${encodeURIComponent(idA)}&b=${encodeURIComponent(idB)}`
          : ''

        cells.push(esc(score), esc(link))
      }

      cells.push(esc(gene.curated_note))
      return cells
    })

    const content = [headers, ...tsvRows].map(r => r.join('\t')).join('\n')
    const blob = new Blob([content], { type: 'text/tab-separated-values' })
    const url = URL.createObjectURL(blob)
    const anchor = document.createElement('a')
    anchor.href = url
    anchor.download = 'transcript_table.tsv'
    anchor.click()
    URL.revokeObjectURL(url)
  }

  return (
    <div className="flex flex-col flex-1 overflow-hidden">
      <div className="bg-white border-b px-4 py-2 flex items-center">
        <button
          onClick={exportTsv}
          className="px-3 py-1.5 text-sm border rounded hover:bg-gray-50 flex items-center gap-1.5"
        >
          ↓ Export TSV ({rows.length} rows)
        </button>
      </div>
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
