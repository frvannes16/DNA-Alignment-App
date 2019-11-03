from django.http import JsonResponse
from django.core import serializers
from dataclasses import dataclass
from .search import DnaSearchTool
from .models import Search, Result
import json


def no_path(request):
    return JsonResponse({'message': 'Welcome to the DNA alignment tool API'})


def submit_search(request):
    if request.method == 'POST':
        try:
            json_request = json.loads(request.body.decode('utf8'))
        except Exception as e:
            return JsonResponse({'submissionMessages': [
                {'type': 'ERROR',
                    'message': 'Invalid search string. JSON was interrupted: ' + str(e)}
            ]})

        if 'searchString' in json_request:
            dna_query = json_request['searchString']
            # TODO: perform validation and string cleaning.
            DnaSearchTool.start_search(dna_query)
        return JsonResponse({'submissionMessages': [
            {'type': 'SUCCESS', 'message': 'Search Started'}
        ]})
    else:
        return JsonResponse({'submissionMessages': [
            {'type': 'ERROR', 'message': 'This endpoint only accepts AJAX POST request.'}
        ]})


def poll_searches(request):
    all_searches = Search.objects.all()

    return JsonResponse({'message': 'endpoint not supported'}, status=404)
