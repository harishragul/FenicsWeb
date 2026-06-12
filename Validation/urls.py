from django.urls import path

from .views import benchmark, convergence, dashboard

urlpatterns = [
    # 'validation' is the pre-existing URL name referenced by every page's
    # navbar (was previously a basic 1D-only check in SteadyStateThermal);
    # it now points at the comprehensive Validation dashboard.
    path('', dashboard, name='validation'),
    path('convergence/', convergence, name='validation_convergence'),
    path('benchmark/', benchmark, name='validation_benchmark'),
]
