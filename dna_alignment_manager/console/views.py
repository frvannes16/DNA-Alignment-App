from django.shortcuts import render
from django.http import HttpResponse

# Create your views here.
def index(request):
    return HttpResponse(b'<html><body><h1>DNA Search Tool</h1></body></html>')

