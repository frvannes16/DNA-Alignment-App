import * as React from 'react';

export enum MessageType {
  ERROR,
  SUCCESS
}

export interface SubmissionMessage {
  type: MessageType;
  message: string;
}

export interface SubmissionResponse {
  submissionMessages: Array<SubmissionMessage>;
}

export interface SubmissionRequest {
  searchString: string;
}
