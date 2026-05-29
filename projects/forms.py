from django import forms
from .models import Project

MESH_TYPE_CHOICES = [
    ('interval',      'IntervalMesh (1D)'),
    ('unit_interval', 'UnitIntervalMesh (1D)'),
    ('rectangle',     'RectangleMesh (2D)'),
    ('unit_square',   'UnitSquareMesh (2D)'),
    ('box',           'BoxMesh (3D)'),
    ('unit_cube',     'UnitCubeMesh (3D)'),
]

BC_TYPES = [
    ('dirichlet', 'Dirichlet (fixed temperature)'),
    ('neumann',   'Neumann (heat flux)'),
    ('robin',     'Robin (convective cooling)'),
]

DIMENSIONS = [('1D', '1D'), ('2D', '2D'), ('3D', '3D')]

_ctrl = {'class': 'form-control'}
_sel  = {'class': 'form-select'}


class ProjectCreateForm(forms.ModelForm):
    class Meta:
        model  = Project
        fields = ['name', 'solver_type']
        widgets = {
            'name':        forms.TextInput(attrs={**_ctrl, 'placeholder': 'e.g. Rod Analysis'}),
            'solver_type': forms.Select(attrs=_sel),
        }


class MeshStepForm(forms.Form):
    mesh_type = forms.ChoiceField(
        choices=MESH_TYPE_CHOICES,
        widget=forms.Select(attrs={**_sel, 'id': 'id_mesh_type'}),
        label='Mesh Type',
    )
    # 1D fields
    n  = forms.IntegerField(required=False, label='Number of cells (n)', min_value=1,
                             widget=forms.NumberInput(attrs=_ctrl))
    x0 = forms.FloatField(required=False, label='x₀ (start)', initial=0.0,
                           widget=forms.NumberInput(attrs=_ctrl))
    x1 = forms.FloatField(required=False, label='x₁ (end)',   initial=1.0,
                           widget=forms.NumberInput(attrs=_ctrl))
    # 2D extra
    y0 = forms.FloatField(required=False, label='y₀ (bottom)', initial=0.0,
                           widget=forms.NumberInput(attrs=_ctrl))
    y1 = forms.FloatField(required=False, label='y₁ (top)',    initial=1.0,
                           widget=forms.NumberInput(attrs=_ctrl))
    nx = forms.IntegerField(required=False, label='nx', min_value=1,
                             widget=forms.NumberInput(attrs=_ctrl))
    ny = forms.IntegerField(required=False, label='ny', min_value=1,
                             widget=forms.NumberInput(attrs=_ctrl))
    # 3D extra
    z0 = forms.FloatField(required=False, label='z₀ (front)', initial=0.0,
                           widget=forms.NumberInput(attrs=_ctrl))
    z1 = forms.FloatField(required=False, label='z₁ (back)',  initial=1.0,
                           widget=forms.NumberInput(attrs=_ctrl))
    nz = forms.IntegerField(required=False, label='nz', min_value=1,
                             widget=forms.NumberInput(attrs=_ctrl))

    def clean(self):
        cleaned = super().clean()
        mt = cleaned.get('mesh_type', '')

        # Required fields per mesh type — only check what the user can actually fill in.
        REQUIRED = {
            'interval':      ['n', 'x0', 'x1'],
            'unit_interval': ['n'],
            'rectangle':     ['nx', 'ny', 'x0', 'x1', 'y0', 'y1'],
            'unit_square':   ['nx', 'ny'],
            'box':           ['nx', 'ny', 'nz', 'x0', 'x1', 'y0', 'y1', 'z0', 'z1'],
            'unit_cube':     ['nx', 'ny', 'nz'],
        }
        for f in REQUIRED.get(mt, []):
            val = cleaned.get(f)
            if val is None and f not in self.errors:
                # Use add_error on the field itself so the input is highlighted red
                label = self.fields[f].label or f
                self.add_error(f, f'This field is required for {mt}.')
        return cleaned


def _bc_fields(prefix):
    """Return a flat dict of BC fields for one face, prefixed."""
    return {
        f'{prefix}_type':  forms.ChoiceField(
            choices=BC_TYPES,
            label=f'{prefix.capitalize()} BC type',
            initial='dirichlet',
            widget=forms.Select(attrs={**_sel, 'class': 'form-select bc-type-select',
                                       'data-face': prefix}),
        ),
        f'{prefix}_value': forms.FloatField(
            required=False, initial=0.0, label='Temperature / Flux value',
            widget=forms.NumberInput(attrs=_ctrl),
        ),
        f'{prefix}_h':     forms.FloatField(
            required=False, label='Convection coefficient h',
            widget=forms.NumberInput(attrs=_ctrl),
        ),
        f'{prefix}_u_inf': forms.FloatField(
            required=False, label='Ambient temperature T∞',
            widget=forms.NumberInput(attrs=_ctrl),
        ),
    }


class SetupStepForm(forms.Form):
    # dimension is intentionally excluded — it is inferred from the mesh type in the view
    k = forms.FloatField(label='Thermal conductivity k (W/m·K)', initial=1.0,
                          widget=forms.NumberInput(attrs=_ctrl))
    f = forms.FloatField(label='Volumetric source term f',        initial=0.0,
                          widget=forms.NumberInput(attrs=_ctrl))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for face in ('left', 'right', 'bottom', 'top', 'front', 'back'):
            self.fields.update(_bc_fields(face))
