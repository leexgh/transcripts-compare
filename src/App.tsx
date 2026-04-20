import React, { useMemo } from 'react'
import { HashRouter, Routes, Route, NavLink, useNavigate, useLocation, useSearchParams } from 'react-router-dom'
import { useGeneData, useFilteredGenes } from './hooks/useGeneData'
import FilterBar from './components/FilterBar'
import TranscriptTable from './components/TranscriptTable'
import DiffViewer from './components/DiffViewer'
import ListCompareTab from './components/ListCompareTab'
import TranscriptCompareTab from './components/TranscriptCompareTab'
import ProteinComparePage from './components/ProteinComparePage'
import type { FilterState, AssemblyFilter } from './lib/types'
import type { SortingState } from '@tanstack/react-table'

function MainView() {
  const { data, loading, error } = useGeneData()
  const [searchParams, setSearchParams] = useSearchParams()

  const filters: FilterState = {
    assembly: (['both', 'GRCh37', 'GRCh38'].includes(searchParams.get('assembly') ?? '')
      ? searchParams.get('assembly') as AssemblyFilter : 'both'),
    search:       searchParams.get('q') ?? '',
    clinicalOnly:  searchParams.has('clinical'),
    germlineOnly:  searchParams.has('germline'),
    maneOnly:      searchParams.has('mane'),
    mismatchOnly:  searchParams.has('mismatch'),
  }

  const sorting: SortingState = useMemo(() => {
    const raw = searchParams.get('sort')
    if (!raw) return []
    return raw.split(',').flatMap(part => {
      const [id, dir] = part.split(':')
      return id ? [{ id, desc: dir === 'desc' }] : []
    })
  }, [searchParams])

  const setFilters = (f: FilterState) => {
    setSearchParams(prev => {
      const p = new URLSearchParams(prev)
      if (f.assembly !== 'both') p.set('assembly', f.assembly); else p.delete('assembly')
      if (f.search)       p.set('q', f.search);    else p.delete('q')
      if (f.clinicalOnly) p.set('clinical', '1');  else p.delete('clinical')
      if (f.germlineOnly) p.set('germline', '1');  else p.delete('germline')
      if (f.maneOnly)     p.set('mane', '1');      else p.delete('mane')
      if (f.mismatchOnly) p.set('mismatch', '1');  else p.delete('mismatch')
      return p
    }, { replace: true })
  }

  const setSorting = (updater: SortingState | ((prev: SortingState) => SortingState)) => {
    const next = typeof updater === 'function' ? updater(sorting) : updater
    setSearchParams(prev => {
      const p = new URLSearchParams(prev)
      if (next.length === 0) p.delete('sort')
      else p.set('sort', next.map(s => `${s.id}:${s.desc ? 'desc' : 'asc'}`).join(','))
      return p
    }, { replace: true })
  }

  const genes    = data?.genes ?? []
  const filtered = useFilteredGenes(genes, filters)

  if (loading) return <div className="flex-1 flex items-center justify-center text-gray-500">Loading data…</div>
  if (error)   return <div className="flex-1 flex items-center justify-center text-red-500">Error: {error}</div>
  if (!data)   return null

  return (
    <div className="flex flex-col flex-1 overflow-hidden">
      <FilterBar filters={filters} onChange={setFilters} geneCount={filtered.length} total={genes.length} />
      <TranscriptTable data={data} genes={filtered} sorting={sorting} onSortingChange={setSorting} />
    </div>
  )
}

function DiffPage() {
  const { data } = useGeneData()
  return (
    <div className="flex-1 overflow-auto">
      <DiffViewer data={data} />
    </div>
  )
}

function ComparePage() {
  const { data, loading, error } = useGeneData()
  if (loading) return <div className="flex-1 flex items-center justify-center text-gray-500">Loading…</div>
  if (error || !data) return <div className="flex-1 p-4 text-red-500">{error ?? 'No data'}</div>
  return <ListCompareTab data={data} />
}

function TranscriptComparePage() {
  const { data, loading, error } = useGeneData()
  if (loading) return <div className="flex-1 flex items-center justify-center text-gray-500">Loading data…</div>
  if (error)   return <div className="flex-1 flex items-center justify-center text-red-500">Error: {error}</div>
  if (!data)   return null
  return <TranscriptCompareTab data={data} />
}

function Nav() {
  const loc = useLocation()
  const tabClass = (active: boolean) =>
    `px-4 py-2 text-sm font-medium border-b-2 ${active ? 'border-blue-600 text-blue-600' : 'border-transparent text-gray-600 hover:text-gray-800'}`

  return (
    <nav className="bg-white border-b px-4 flex items-center gap-0">
      <span className="font-bold text-gray-800 mr-6">Transcript Diff</span>
      <NavLink to="/table" className={({ isActive }) => tabClass(isActive)}>
        Transcript Table
      </NavLink>
      <NavLink to="/" end className={({ isActive }) => tabClass(isActive)}>
        Multiple Transcripts Compare
      </NavLink>
      <NavLink to="/compare" className={({ isActive }) => tabClass(isActive)}>
        List Compare
      </NavLink>
    </nav>
  )
}

export default function App() {
  return (
    <HashRouter>
      <div className="h-screen flex flex-col bg-gray-50">
        <Routes>
          <Route path="/diff" element={
            <>
              <Nav />
              <DiffPage />
            </>
          } />
          <Route path="/compare-protein" element={
            <ProteinComparePage />
          } />
          <Route path="/compare" element={
            <>
              <Nav />
              <ComparePage />
            </>
          } />
          <Route path="/table" element={
            <>
              <Nav />
              <MainView />
            </>
          } />
          <Route path="*" element={
            <>
              <Nav />
              <TranscriptComparePage />
            </>
          } />
        </Routes>
      </div>
    </HashRouter>
  )
}
