import * as React from 'react';
import axios from 'axios';
import { Result } from './types';
import Results from './Results';

const POLL_INTERVAL = 3000; // ms

/** The ResultContainer polls for search results every 3 seconds
 * and provides the results to Results.tsx for presentation.
 */
const ResultsContainer = () => {
  const [results, setResults] = React.useState<Array<Result>>();

  const pollResults = () => {
    axios.get<Array<Result>>('/api/search/poll').then(response => {
      if (response && response.data && response.data) {
        console.log('Poll response recevied');
        setResults(response.data);
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
