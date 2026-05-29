from django.urls import path
from .views import conduction, validation

urlpatterns = [
    path('', conduction, name='conduction'),
    path('validation/', validation, name='validation'),
]
