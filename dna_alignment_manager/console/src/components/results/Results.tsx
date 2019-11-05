import * as React from 'react';
import { Result } from './types';

interface Props {
  results: Array<Result>;
}

const Results = (props: Props) =>
  <table>
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
            <tr key={result.searchId}>
              <td>
                {result.searchId}
              </td>
              <td>
                {result.searchString}
              </td>
              <td>
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
