import React, { useState } from 'react'
import { HashRouter, Routes, Route, NavLink, useNavigate, useLocation } from 'react-router-dom'
import { useGeneData, useFilteredGenes } from './hooks/useGeneData'
import FilterBar from './components/FilterBar'
import TranscriptTable from './components/TranscriptTable'
import DiffViewer from './components/DiffViewer'
import ListCompareTab from './components/ListCompareTab'
import TranscriptCompareTab from './components/TranscriptCompareTab'
import ProteinComparePage from './components/ProteinComparePage'
import type { FilterState } from './lib/types'

const DEFAULT_FILTERS: FilterState = {
  assembly:     'both',
  search:       '',
  clinicalOnly:  false,
  germlineOnly:  false,
  maneOnly:      false,
  mismatchOnly:  false,
}

function MainView() {
  const { data, loading, error } = useGeneData()
  const [filters, setFilters]   = useState<FilterState>(DEFAULT_FILTERS)

  const genes    = data?.genes ?? []
  const filtered = useFilteredGenes(genes, filters)

  if (loading) return <div className="flex-1 flex items-center justify-center text-gray-500">Loading data…</div>
  if (error)   return <div className="flex-1 flex items-center justify-center text-red-500">Error: {error}</div>
  if (!data)   return null

  return (
    <div className="flex flex-col flex-1 overflow-hidden">
      <FilterBar filters={filters} onChange={setFilters} geneCount={filtered.length} total={genes.length} />
      <TranscriptTable data={data} genes={filtered} />
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
