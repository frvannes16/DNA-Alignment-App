import * as React from 'react';
import axios from 'axios';

/** This tells axios where to find the csrf token and to add it as a header. */
axios.defaults.xsrfCookieName = 'csrftoken'
axios.defaults.xsrfHeaderName = "X-CSRFTOKEN"
axios.defaults.withCredentials = true

enum MessageType {
  ERROR,
  SUCCESS
}

interface SubmissionMessage {
  type: MessageType;
  message: string;
}

interface SubmissionResponse {
  submissionMessages: Array<SubmissionMessage>;
}

interface SubmissionRequest {
  searchString: string;
}

/**
 * The Search Tool is responsible for submitting new searches and displaying the status and results of previous search results.
 */
const SearchTool = () => {
  const [searchDna, setSearchDna] = React.useState<string>('');
  const [submissionMessages, setSubmissionMessages] = React.useState<
    Array<SubmissionMessage> | undefined
  >(undefined);

  const executeSearch = () => {
    setSubmissionMessages(null);
    // submit the search to the API.
    if (searchDna.length > 0) {
      const submission: SubmissionRequest = { searchString: searchDna };
      axios.post<SubmissionResponse>('/api/search/submit', submission).then(response => {
        if (response && response.data && response.data.submissionMessages) {
          setSubmissionMessages(response.data.submissionMessages);
          if (!response.data.submissionMessages.some(msg => msg.type === MessageType.ERROR)) {
            setSearchDna('')
          }
        }
      });
    }
  };

  return (
    <div>
      <h3>Submit a new search</h3>
      <div>
        {submissionMessages &&
          submissionMessages.map(message =>
            <div key={message.message}>
              <span>{message.type}:</span><span> {message.message}</span>
            </div>
          )}
      </div>
      <form name="search-form">
        <input
          type="text"
          value={searchDna}
          onChange={e => setSearchDna(e.target.value)}
        />
        <button type="button" onClick={executeSearch}>
          Search
        </button>
      </form>
    </div>
  );
};

export default SearchTool;
