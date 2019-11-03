from django.urls import path

from . import views

urlpatterns = [
    path('', views.no_path, name='no_path'),
    path('search/submit', views.submit_search, name='submit-search'),
    path('search/poll/', views.poll_searches, name='poll-searches')

]