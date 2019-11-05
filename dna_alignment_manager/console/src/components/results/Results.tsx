import * as React from 'react';
import { Result } from './types';

interface Props {
  results: Array<Result>;
}

const getStatusStyle = (status: string) => {
  switch (status) {
    case 'PROCESSING':
      return 'processing';
    case 'FAILED':
      return 'failed';
    case 'MATCH_FOUND':
      return 'matchFound';
    case 'NO_MATCH':
    default:
      return 'noMatch';
  }
};

const Results = (props: Props) =>
  <table className={'results'}>
    <thead>
      <tr>
        <th>ID</th>
        <th>Search String</th>
        <th>Status</th>
        <th>Matched Protein</th>
        <th>Alignment start index</th>
        <th>Alignment end index</th>
      </tr>
    </thead>
    <tbody>
      {props.results &&
        props.results.map(result => {
          return (
            <tr key={result.searchId} className={'row'}>
              <td>
                {result.searchId}
              </td>
              <td>
                {result.searchString}
              </td>
              <td className={`${getStatusStyle(result.status)} status`}>
                {result.status}
              </td>
              <td>
                {result.match && result.match.protein}
              </td>
              <td>
                {result.match && result.match.startPos}
              </td>
              <td>
                {result.match && result.match.endPos}
              </td>
            </tr>
          );
        })}
    </tbody>
  </table>;

export default Results;
