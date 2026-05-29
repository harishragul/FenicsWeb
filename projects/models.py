from django.db import models


class Project(models.Model):
    SOLVER_CHOICES = [
        ('steady_conduction',    'Steady Conduction'),
        ('convection_diffusion', 'Convection-Diffusion'),
        ('transient',            'Transient Heat'),
        ('inverse',              'Inverse Problem'),
    ]
    STATUS = [
        ('new',        'New'),
        ('mesh_done',  'Mesh Complete'),
        ('setup_done', 'Setup Complete'),
        ('solved',     'Solved'),
    ]
    name        = models.CharField(max_length=200)
    solver_type = models.CharField(max_length=50, choices=SOLVER_CHOICES,
                                    default='steady_conduction')
    status      = models.CharField(max_length=20, choices=STATUS, default='new')
    created_at  = models.DateTimeField(auto_now_add=True)
    updated_at  = models.DateTimeField(auto_now=True)

    def __str__(self):
        return f"{self.name} ({self.get_solver_type_display()})"

    @property
    def mesh_plot_url(self):
        if hasattr(self, 'mesh_config') and self.mesh_config.mesh_plot_path:
            return f'/media/{self.mesh_config.mesh_plot_path}'
        return None

    @property
    def solution_plot_url(self):
        if hasattr(self, 'result') and self.result.solution_plot_path:
            return f'/media/{self.result.solution_plot_path}'
        return None

    @property
    def step_status(self):
        return {
            'mesh':  self.status in ('mesh_done', 'setup_done', 'solved'),
            'setup': self.status in ('setup_done', 'solved'),
            'solve': self.status == 'solved',
        }


class ProjectMesh(models.Model):
    project        = models.OneToOneField(Project, on_delete=models.CASCADE,
                                           related_name='mesh_config')
    mesh_type      = models.CharField(max_length=50)
    mesh_str       = models.TextField()
    mesh_params    = models.JSONField()
    mesh_plot_path = models.CharField(max_length=300, blank=True)


class ProjectSetup(models.Model):
    project       = models.OneToOneField(Project, on_delete=models.CASCADE,
                                          related_name='setup_config')
    dimension     = models.CharField(max_length=5)
    bc_config     = models.JSONField()
    solver_params = models.JSONField()


class ProjectResult(models.Model):
    project             = models.OneToOneField(Project, on_delete=models.CASCADE,
                                               related_name='result')
    mesh_plot_path      = models.CharField(max_length=300, blank=True)
    solution_plot_path  = models.CharField(max_length=300, blank=True)
    solution_data       = models.JSONField(default=list)
    solved_at           = models.DateTimeField(auto_now_add=True)
    solve_time_s        = models.FloatField(null=True)
