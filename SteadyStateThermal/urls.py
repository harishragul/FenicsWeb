from django.urls import path
from .views import conduction

urlpatterns = [
    path('', conduction, name='conduction'),
]
