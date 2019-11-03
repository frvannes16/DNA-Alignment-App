import * as React from 'react';

enum MessageType {
  ERROR,
  SUCCESS
}

interface SubmissionMessage {
  type: MessageType;
  message: string;
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
    // If it succeeds, clear the input and poll for results.
    // If it fails, show a failure message and keep the input.
  };

  return (
    <div>
      <h3>Submit a new search</h3>
      <div>
        {submissionMessages &&
          submissionMessages.map(message =>
            <div>
              {message}
            </div>
          )}
      </div>
      <form name="search-form">
        <input
          type="text"
          value={searchDna}
          onChange={e => setSearchDna(e.target.value)}
        />
        <button type="submit" onClick={executeSearch}>
          Search
        </button>
      </form>
    </div>
  );
};

export default SearchTool;
