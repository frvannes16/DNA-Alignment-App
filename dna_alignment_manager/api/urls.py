from django.urls import path

from . import views

urlpatterns = [
    path('', views.no_path, name='no_path'),
]