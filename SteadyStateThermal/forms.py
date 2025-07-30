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
    
    BC_TYPES = [
        ('dirichlet', 'Dirichlet (Fixed Temperature)'),
        ('convection', 'Convection'),
    ]
    
    dimension = forms.ChoiceField(
        choices=DIMENSIONS, 
        label="Select Dimension",
        widget=forms.Select(attrs={'class': 'form-select'})
    )
    
    # Boundary condition types
    left_bc_type = forms.ChoiceField(
        choices=BC_TYPES, 
        label="Left BC Type", 
        initial='dirichlet',
        widget=forms.Select(attrs={'class': 'form-select'})
    )
    right_bc_type = forms.ChoiceField(
        choices=BC_TYPES, 
        label="Right BC Type", 
        initial='dirichlet',
        widget=forms.Select(attrs={'class': 'form-select'})
    )
    top_bc_type = forms.ChoiceField(
        choices=BC_TYPES, 
        label="Top BC Type (for 2D/3D)", 
        initial='dirichlet', 
        required=False,
        widget=forms.Select(attrs={'class': 'form-select'})
    )
    bottom_bc_type = forms.ChoiceField(
        choices=BC_TYPES, 
        label="Bottom BC Type (for 2D/3D)", 
        initial='dirichlet', 
        required=False,
        widget=forms.Select(attrs={'class': 'form-select'})
    )
    front_bc_type = forms.ChoiceField(
        choices=BC_TYPES, 
        label="Front BC Type (for 3D)", 
        initial='dirichlet', 
        required=False,
        widget=forms.Select(attrs={'class': 'form-select'})
    )
    back_bc_type = forms.ChoiceField(
        choices=BC_TYPES, 
        label="Back BC Type (for 3D)", 
        initial='dirichlet', 
        required=False,
        widget=forms.Select(attrs={'class': 'form-select'})
    )
    
    # Dirichlet boundary condition values
    left_bc = forms.FloatField(
        label="Left Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    right_bc = forms.FloatField(
        label="Right Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    top_bc = forms.FloatField(
        label="Top Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    bottom_bc = forms.FloatField(
        label="Bottom Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    front_bc = forms.FloatField(
        label="Front Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    back_bc = forms.FloatField(
        label="Back Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    
    # Convection parameters
    left_h = forms.FloatField(
        label="Left Convection Coefficient (W/m²K)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    right_h = forms.FloatField(
        label="Right Convection Coefficient (W/m²K)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    top_h = forms.FloatField(
        label="Top Convection Coefficient (W/m²K)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    bottom_h = forms.FloatField(
        label="Bottom Convection Coefficient (W/m²K)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    front_h = forms.FloatField(
        label="Front Convection Coefficient (W/m²K)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    back_h = forms.FloatField(
        label="Back Convection Coefficient (W/m²K)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    
    left_t_amb = forms.FloatField(
        label="Left Ambient Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    right_t_amb = forms.FloatField(
        label="Right Ambient Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    top_t_amb = forms.FloatField(
        label="Top Ambient Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    bottom_t_amb = forms.FloatField(
        label="Bottom Ambient Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    front_t_amb = forms.FloatField(
        label="Front Ambient Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    back_t_amb = forms.FloatField(
        label="Back Ambient Temperature (°C)", 
        required=False,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    
    f = forms.FloatField(
        label="Source Term (f)", 
        initial=0.0,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
    k = forms.FloatField(
        label="Thermal Conductivity (k)", 
        initial=1.0,
        widget=forms.NumberInput(attrs={'class': 'form-control', 'step': '0.1'})
    )
