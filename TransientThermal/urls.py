from django.urls import path
from .views import transient, scheme_comparison

urlpatterns = [
    path('', transient, name='transient'),
    path('comparison/', scheme_comparison, name='scheme_comparison'),
]
