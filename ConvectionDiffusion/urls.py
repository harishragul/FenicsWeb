from django.urls import path
from .views import convection, peclet_study

urlpatterns = [
    path('', convection, name='convection'),
    path('peclet/', peclet_study, name='peclet_study'),
]
