# Generated by Django 4.2.11 on 2024-12-05 16:32

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('mesh', '0002_mesh_name'),
    ]

    operations = [
        migrations.AlterField(
            model_name='mesh',
            name='name',
            field=models.CharField(max_length=100),
        ),
    ]
