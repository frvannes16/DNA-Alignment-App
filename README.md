# DNA-Alignment-App

The goal of this project was to produce a complete web app that, given a DNA string, searches a small set of known proteins.
This DNA  alignment app is set up to easily scale horizontally for fast asynchronous procesing of alignments. Unfortunately, the app can only perform a simple string search, as I was unable to figure out how to use the BioPython api to get the HSP without querying a remote or local BLAST+ service within the given time constraint.

## How it works
You can view the app at dev.franklinvannes.com:3000/ - This is my home dev server, so please be aware that I might have deactivated the service and closed those ports at the time of reading.

The front end is written in React and Typescript, and it uses two endpoints to communicate with the backend:

### /api/search/submit
New searches are submitted through this endpoint. The search request is stored in the database and a set of search tasks are issued to a distributed and scalable set of celery workers with a Redis broker. Each celery task that is provided to a worker is to search the provided sequence against a protein sequence. If I had the DNA alignment working, then this asynchronous processing by worker threads would be an excellent way to distribute processing of search tasks. Once a worker finishes processing the match details are saved to the database and are compared against the other results. Once the result with the highest score has been determined, that result is marked as the best and chosen result. 

### /api/search/poll
This endpoint is queried every second, and it simply queries the database for searches, and checks if all of the results have been returned. If they have, the endpoint reports the cached best result to the user. If no result has been found or if the results are still processing, the front end is notified in the poll response.

## MVP
This is an MVP that is running on my home server. It would quickly scale onto a computing cluster, and I began a `deployment` branch that contains kubernetes manifest files that would allow for a k8s deployment with a little more tweaking.

