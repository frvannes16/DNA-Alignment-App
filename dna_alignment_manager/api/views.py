from django.http import JsonResponse
from django.core import serializers
from dataclasses import dataclass
import json


def no_path(request):
    return JsonResponse({'message': 'Welcome to the DNA alignment tool API'})


def submit_search(request):
    if request.is_ajax() and request.method == 'POST':
        try:
            json_request = json.load(request.body)
        except Exception as e:
            return JsonResponse({'submissionMessages': [
                {'type': 'ERROR',
                    'message': 'Invalid search string. JSON was interrupted: ' + str(e)}
            ]})

        if 'searchString' in json_request:
            dna_query = json_request['searchString']
            # TODO: perform validation and string cleaning.
            # Start search process. Ensure process has started.
        return JsonResponse({'submissionMessages': [
            {'type': 'SUCCESS', 'message': 'Search Started'}
        ]})
    else:
        return JsonResponse({'submissionMessages': [
            {'type': 'ERROR', 'message': 'This endpoint only accepts AJAX POST request.'}
        ]})


def poll_searches(request):
    return JsonResponse({'message': 'endpoint not supported'}, status=404)
