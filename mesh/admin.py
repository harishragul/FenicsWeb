from django.contrib import admin
from mesh.models import Mesh
# Register your models here.

class MeshAdmin(admin.ModelAdmin):
    list_display = ["name", "mesh", "id"]

admin.site.register(Mesh, MeshAdmin)