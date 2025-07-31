from django.urls import path
from .views import *

urlpatterns = [
    path('', conduction, name='conduction'),
    path('download-report/', generate_pdf_report, name='download_pdf_report'),
]
