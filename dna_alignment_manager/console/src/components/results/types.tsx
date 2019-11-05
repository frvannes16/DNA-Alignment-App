export interface Match {
    protein: string;
    startPos: string;
    endPos: string;
}

export interface Result {
    searchId: number;
    status: string;
    searchString: string;
    match?: Match;    
}

export interface SearchPollResponse {
    searches: Array<Result>;
}