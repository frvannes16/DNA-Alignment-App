import * as React from 'react';
import { Result } from './types';

interface Props {
  results: Array<Result>;
}

const Results = (props: Props) =>
  <div>
    {props.results.map(result => {
      return (
        <div>
          {result.status}
        </div>
      );
    })}
  </div>;

export default Results;
