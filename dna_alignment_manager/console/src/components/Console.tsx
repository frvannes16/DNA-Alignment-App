import * as React from 'react'
import SearchTool from './search-tool/SearchTool';
import ResultsContainer from './results/ResultsContainer';

const Console = () =>
  <div className={'content'}>
    <h1>DNA Search Tool</h1>
    <h2>Franklin van Nes</h2>
    <SearchTool/>
    <h2>Results</h2>
    <ResultsContainer />
  </div>

export default Console
