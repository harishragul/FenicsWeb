from django import forms
from mesh.models import Mesh

BC_TYPES = [
    ('dirichlet', 'Dirichlet (fixed temperature)'),
    ('neumann', 'Neumann (heat flux)'),
    ('robin', 'Robin (convective cooling)'),
]

DIMENSIONS = [
    ('1D', '1D'),
    ('2D', '2D'),
    ('3D', '3D'),
]


class MeshForm(forms.Form):
    mesh = forms.ModelChoiceField(
        queryset=Mesh.objects.all(),
        label="Select the Mesh",
        widget=forms.Select(attrs={'class': 'form-select'}),
    )


class HeatSolverForm(forms.Form):
    dimension = forms.ChoiceField(choices=DIMENSIONS, label="Dimension")

    # Left face
    left_bc_type = forms.ChoiceField(
        choices=BC_TYPES, label="Left BC Type", initial='dirichlet',
        widget=forms.Select(attrs={'class': 'form-select bc-type-select', 'data-face': 'left'}),
    )
    left_bc = forms.FloatField(label="Left: Temperature or Flux value", required=False, initial=0.0)
    left_h = forms.FloatField(label="Left: Convection coefficient h (Robin)", required=False)
    left_u_inf = forms.FloatField(label="Left: Ambient temperature T∞ (Robin)", required=False)

    # Right face
    right_bc_type = forms.ChoiceField(
        choices=BC_TYPES, label="Right BC Type", initial='dirichlet',
        widget=forms.Select(attrs={'class': 'form-select bc-type-select', 'data-face': 'right'}),
    )
    right_bc = forms.FloatField(label="Right: Temperature or Flux value", required=False, initial=0.0)
    right_h = forms.FloatField(label="Right: Convection coefficient h (Robin)", required=False)
    right_u_inf = forms.FloatField(label="Right: Ambient temperature T∞ (Robin)", required=False)

    # Bottom face (2D / 3D)
    bottom_bc_type = forms.ChoiceField(
        choices=BC_TYPES, label="Bottom BC Type", initial='dirichlet',
        widget=forms.Select(attrs={'class': 'form-select bc-type-select', 'data-face': 'bottom'}),
    )
    bottom_bc = forms.FloatField(label="Bottom: Temperature or Flux value", required=False, initial=0.0)
    bottom_h = forms.FloatField(label="Bottom: h (Robin)", required=False)
    bottom_u_inf = forms.FloatField(label="Bottom: T∞ (Robin)", required=False)

    # Top face (2D / 3D)
    top_bc_type = forms.ChoiceField(
        choices=BC_TYPES, label="Top BC Type", initial='dirichlet',
        widget=forms.Select(attrs={'class': 'form-select bc-type-select', 'data-face': 'top'}),
    )
    top_bc = forms.FloatField(label="Top: Temperature or Flux value", required=False, initial=0.0)
    top_h = forms.FloatField(label="Top: h (Robin)", required=False)
    top_u_inf = forms.FloatField(label="Top: T∞ (Robin)", required=False)

    # Front face (3D only)
    front_bc_type = forms.ChoiceField(
        choices=BC_TYPES, label="Front BC Type", initial='dirichlet',
        widget=forms.Select(attrs={'class': 'form-select bc-type-select', 'data-face': 'front'}),
    )
    front_bc = forms.FloatField(label="Front: Temperature or Flux value", required=False, initial=0.0)
    front_h = forms.FloatField(label="Front: h (Robin)", required=False)
    front_u_inf = forms.FloatField(label="Front: T∞ (Robin)", required=False)

    # Back face (3D only)
    back_bc_type = forms.ChoiceField(
        choices=BC_TYPES, label="Back BC Type", initial='dirichlet',
        widget=forms.Select(attrs={'class': 'form-select bc-type-select', 'data-face': 'back'}),
    )
    back_bc = forms.FloatField(label="Back: Temperature or Flux value", required=False, initial=0.0)
    back_h = forms.FloatField(label="Back: h (Robin)", required=False)
    back_u_inf = forms.FloatField(label="Back: T∞ (Robin)", required=False)

    # Solver parameters
    f = forms.FloatField(label="Source Term f (volumetric heat generation)", initial=0.0)
    k = forms.FloatField(label="Thermal Conductivity k (W/m·K)", initial=1.0)
