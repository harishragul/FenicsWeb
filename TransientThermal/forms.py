from django import forms

DIMENSIONS = [
    ('1D', '1D — Unit Interval'),
    ('2D', '2D — Unit Square'),
]

SCHEMES = [
    ('both', 'Compare Backward Euler vs Crank-Nicolson'),
    ('backward_euler', 'Backward Euler only'),
    ('crank_nicolson', 'Crank-Nicolson only'),
]

INITIAL_CONDITIONS = [
    ('sine', 'Sine mode — sin(πx) [analytical benchmark]'),
    ('gaussian', 'Gaussian pulse centered at 0.5'),
    ('constant', 'Uniform initial temperature (= 1)'),
]

MAX_TIME_STEPS = 200


class TransientForm(forms.Form):
    dimension = forms.ChoiceField(choices=DIMENSIONS, label="Dimension", initial='1D')
    n = forms.IntegerField(
        label="Mesh resolution (cells per side)", initial=20, min_value=2, max_value=50,
    )

    rho = forms.FloatField(label="Density ρ", initial=1.0)
    cp = forms.FloatField(label="Specific heat capacity cp", initial=1.0)
    k = forms.FloatField(label="Conductivity k", initial=1.0)
    f = forms.FloatField(label="Source term f", initial=0.0)

    T = forms.FloatField(label="End time T", initial=0.1)
    dt = forms.FloatField(label="Time step dt", initial=0.01)
    initial_condition = forms.ChoiceField(choices=INITIAL_CONDITIONS, label="Initial condition", initial='sine')
    scheme = forms.ChoiceField(choices=SCHEMES, label="Time integration scheme", initial='both')

    left_bc = forms.FloatField(label="Left boundary value", initial=0.0)
    right_bc = forms.FloatField(label="Right boundary value", initial=0.0)
    bottom_bc = forms.FloatField(label="Bottom boundary value (2D only)", initial=0.0, required=False)
    top_bc = forms.FloatField(label="Top boundary value (2D only)", initial=0.0, required=False)

    def clean(self):
        cleaned_data = super().clean()
        for field in ('bottom_bc', 'top_bc'):
            if cleaned_data.get(field) is None:
                cleaned_data[field] = 0.0

        T = cleaned_data.get('T')
        dt = cleaned_data.get('dt')
        if T is not None and dt is not None:
            if dt <= 0 or T <= 0:
                raise forms.ValidationError("T and dt must be positive.")
            if dt > T:
                raise forms.ValidationError("Time step dt must not exceed end time T.")
            if round(T / dt) > MAX_TIME_STEPS:
                raise forms.ValidationError(
                    f"T/dt = {T / dt:.0f} exceeds the maximum of {MAX_TIME_STEPS} time steps."
                )
        return cleaned_data
