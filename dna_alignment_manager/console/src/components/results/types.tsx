export interface Match {
    protein: string;
    start_pos: string;
    end_pos: string;
    match_confidence: number;
}

export interface Result {
    searchId: number;
    status: string;
    search_string: string;
    match?: Match;    
}