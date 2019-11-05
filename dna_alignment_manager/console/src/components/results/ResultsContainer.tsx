import * as React from 'react';
import axios from 'axios';
import { Result, SearchPollResponse } from './types';
import Results from './Results';

const POLL_INTERVAL = 1000; // ms

/** The ResultContainer polls for search results every POLL_INTERVAL seconds
 * and provides the results to Results.tsx for presentation.
 */
const ResultsContainer = () => {
  const [results, setResults] = React.useState<Array<Result>>([]);

  const pollResults = () => {
    axios.get<SearchPollResponse>('/api/search/poll').then(response => {
      if (response && response.data && response.data.searches) {
        console.log('Poll response received');
        setResults(response.data.searches);
      }
    });
  };

  // Start the polling interval
  React.useEffect(() => {
    const intervalId = setInterval(pollResults, POLL_INTERVAL);
    return () => {
      clearInterval(intervalId);
    };
  });

  return <Results results={results} />;
};

export default ResultsContainer;
