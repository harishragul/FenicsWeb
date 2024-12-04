from django.urls import path
from .views import *

urlpatterns = [
    path('', generate_mesh, name='generate_mesh'),
    path('show_mesh/', show_mesh, name='show_mesh'),
]
