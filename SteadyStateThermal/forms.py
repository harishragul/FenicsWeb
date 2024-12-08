from django import forms
from mesh.models import Mesh

class MeshForm(forms.Form):
    mesh = forms.ModelChoiceField(queryset=Mesh.objects.all(),label="Select the Mesh",widget=forms.Select(attrs={'class': 'form-select'}))

class HeatSolverForm(forms.Form):
    DIMENSIONS = [
        ('1D', '1D'),
        ('2D', '2D'),
        ('3D', '3D'),
    ]
    
    dimension = forms.ChoiceField(choices=DIMENSIONS, label="Select Dimension")
    left_bc = forms.FloatField(label="Left Boundary Condition", required=False)
    right_bc = forms.FloatField(label="Right Boundary Condition", required=False)
    top_bc = forms.FloatField(label="Top Boundary Condition (for 2D/3D)", required=False)
    bottom_bc = forms.FloatField(label="Bottom Boundary Condition (for 2D/3D)", required=False)
    front_bc = forms.FloatField(label="Front Boundary Condition (for 3D)", required=False)
    back_bc = forms.FloatField(label="Back Boundary Condition (for 3D)", required=False)
    f = forms.FloatField(label="Source Term (f)", initial=0.0)
    k = forms.FloatField(label="Thermal Conductivity (k)", initial=1.0)
