from django.shortcuts import render
from django.http import JsonResponse


def no_path(request):
    return JsonResponse({'message': 'Welcome to the DNA alignment tool API'})
