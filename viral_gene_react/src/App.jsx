import React,{useEffect,useState} from 'react'
import { HashRouter as Router, Routes, Route } from 'react-router-dom'
import Home from './components/Home'
import Result from './components/Result'
const App = () => {
  return (
    <>
    <Router>
      <Routes>
        <Route path="/" element={<Home />} />
        <Route path='/result' element={<Result />} />
      </Routes>
    </Router>
    </>
  )
}

export default App
