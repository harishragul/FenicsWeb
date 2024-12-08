from django.urls import path
from .views import *

urlpatterns = [
    path('', conduction, name='conduction'),
]
