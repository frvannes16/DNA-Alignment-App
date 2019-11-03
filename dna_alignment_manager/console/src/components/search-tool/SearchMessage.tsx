import * as React from 'react';
import { SubmissionMessage } from './types';

interface Props {
  messages: Array<SubmissionMessage>;
}

const SearchMessages = (props: Props) =>
  <div>
    {props.messages &&
      props.messages.map(message =>
        <div key={message.message}>
          <span>
            {message.type}:
          </span>
          <span>
            {' '}{message.message}
          </span>
        </div>
      )}
  </div>;
export default SearchMessages;
