from django import forms

DIMENSIONS = [
    ('1D', '1D — Unit Interval'),
    ('2D', '2D — Unit Square'),
]

METHODS = [
    ('both', 'Compare Galerkin vs SUPG'),
    ('galerkin', 'Standard Galerkin only'),
    ('supg', 'SUPG stabilized only'),
]


class ConvectionDiffusionForm(forms.Form):
    dimension = forms.ChoiceField(choices=DIMENSIONS, label="Dimension", initial='1D')
    n = forms.IntegerField(
        label="Mesh resolution (cells per side)", initial=20, min_value=2, max_value=200,
    )

    k = forms.FloatField(label="Diffusivity k", initial=0.01)
    bx = forms.FloatField(label="Velocity bx", initial=1.0)
    by = forms.FloatField(label="Velocity by (2D only)", initial=0.0, required=False)
    f = forms.FloatField(label="Source term f", initial=0.0)

    left_bc = forms.FloatField(label="Left boundary value", initial=0.0)
    right_bc = forms.FloatField(label="Right boundary value", initial=1.0)
    bottom_bc = forms.FloatField(label="Bottom boundary value (2D only)", initial=0.0, required=False)
    top_bc = forms.FloatField(label="Top boundary value (2D only)", initial=0.0, required=False)

    method = forms.ChoiceField(choices=METHODS, label="Method", initial='both')

    def clean(self):
        cleaned_data = super().clean()
        # Default the optional 2D-only fields when left blank.
        for field in ('by', 'bottom_bc', 'top_bc'):
            if cleaned_data.get(field) is None:
                cleaned_data[field] = 0.0
        return cleaned_data
