from django.db import models

class Mesh(models.Model):
    name = models.CharField(max_length=100, null=False, blank=False)
    mesh = models.CharField(max_length=200, null=False, blank=False)

    def __str__(self):
        return f"{self.name} - {self.id}"