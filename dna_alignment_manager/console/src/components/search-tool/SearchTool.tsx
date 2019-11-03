import * as React from 'react';

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

const apiSearch = (
  submission: SubmissionRequest
): Promise<SubmissionResponse> => {
  const url = '/api/search/submit';
  return fetch(url, {
    method: 'POST',
    body: JSON.stringify(submission),
    cache: 'no-cache',
    headers: {
      'Content-Type': 'application/json'
    }
  }).then(response => {
    if (!response.ok) {
      throw new Error(response.statusText);
    }
    return response.json() as Promise<SubmissionResponse>;
  });
};

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
      apiSearch(submission).then(response => {
        setSubmissionMessages(response.submissionMessages);
        if (
          !response.submissionMessages.some(
            msg => msg.type === MessageType.ERROR
          )
        ) {
          // Clear the DNA search field if there were no errors.
          setSearchDna('');
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
