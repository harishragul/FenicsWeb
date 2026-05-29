from django.urls import path
from . import views

urlpatterns = [
    path('',                 views.home,       name='home'),
    path('new/',             views.create,     name='project_create'),
    path('<int:pk>/',        views.overview,   name='project_overview'),
    path('<int:pk>/mesh/',   views.mesh_step,  name='project_mesh'),
    path('<int:pk>/setup/',  views.setup_step, name='project_setup'),
    path('<int:pk>/solve/',  views.solve_step, name='project_solve'),
    path('<int:pk>/delete/', views.delete,     name='project_delete'),
]
