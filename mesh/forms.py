from django import forms

# Choices for mesh types
MESH_TYPE_CHOICES = [
    ('interval', 'IntervalMesh (1D)'),
    ('rectangle', 'RectangleMesh (2D)'),
    ('box', 'BoxMesh (3D)'),
    ('unit_interval', 'UnitIntervalMesh (1D)'),
    ('unit_square', 'UnitSquareMesh (2D)'),
    ('unit_cube', 'UnitCubeMesh (3D)'),
    ('circle', 'CircleMesh (2D)'),
    ('sphere', 'SphereMesh (3D)'),
    #('custom', 'Custom Mesh (Gmsh, XML, XDMF)'),
]

# MeshType
class MeshType(forms.Form):
    mesh_name = forms.CharField(max_length=200, required=True, label="Mesh Name")
    mesh_type = forms.ChoiceField(choices=MESH_TYPE_CHOICES, label="Mesh Type", required=True)

# IntervalMesh fields
class IntervalMesh(forms.Form):
    interval_n = forms.IntegerField(label="Number of intervals (n)", required=True)
    interval_x0 = forms.FloatField(label="Start point (x0)", required=True)
    interval_x1 = forms.FloatField(label="End point (x1)", required=True)

# RectangleMesh fields
class RectangleMesh(forms.Form):
    rectangle_x0 = forms.FloatField(label="Bottom-left corner x (x0)", required=True)
    rectangle_y0 = forms.FloatField(label="Bottom-left corner y (y0)", required=True)
    rectangle_x1 = forms.FloatField(label="Top-right corner x (x1)", required=True)
    rectangle_y1 = forms.FloatField(label="Top-right corner y (y1)", required=True)
    rectangle_nx = forms.IntegerField(label="Number of divisions in x (nx)", required=True)
    rectangle_ny = forms.IntegerField(label="Number of divisions in y (ny)", required=True)

# BoxMesh fields
class BoxMesh(forms.Form):
    box_x0 = forms.FloatField(label="Bottom-left corner x (x0)", required=True)
    box_y0 = forms.FloatField(label="Bottom-left corner y (y0)", required=True)
    box_z0 = forms.FloatField(label="Bottom-left corner z (z0)", required=True)
    box_x1 = forms.FloatField(label="Top-right corner x (x1)", required=True)
    box_y1 = forms.FloatField(label="Top-right corner y (y1)", required=True)
    box_z1 = forms.FloatField(label="Top-right corner z (z1)", required=True)
    box_nx = forms.IntegerField(label="Number of divisions in x (nx)", required=True)
    box_ny = forms.IntegerField(label="Number of divisions in y (ny)", required=True)
    box_nz = forms.IntegerField(label="Number of divisions in z (nz)", required=True)

# UnitMesh fields
class UnitIntervalMesh(forms.Form):
    unit_nx = forms.IntegerField(label="Number of divisions (nx)", required=True)

class UnitSquareMesh(forms.Form):
    unit_nx = forms.IntegerField(label="Number of divisions (nx)", required=True)
    unit_ny = forms.IntegerField(label="Number of divisions (ny)", required=True)

class UnitCubeMesh(forms.Form):
    unit_nx = forms.IntegerField(label="Number of divisions (nx)", required=True)
    unit_ny = forms.IntegerField(label="Number of divisions (ny)", required=True)
    unit_nz = forms.IntegerField(label="Number of divisions (nz)", required=True)

# ExternalMesh fields
"""
class CustomMesh(forms.Form):
    mesh = forms.FileField(required=True)

class DynamicMeshForm(forms.Form):
    mesh_type = forms.ChoiceField(choices=MESH_TYPE_CHOICES, label="Mesh Type", required=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['mesh_details'] = forms.JSONField(widget=forms.HiddenInput(), required=False)
"""