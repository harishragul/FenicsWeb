from django.urls import path
from .views import *

urlpatterns = [
    path('', generate_mesh_view, name='generate_mesh'),
]
