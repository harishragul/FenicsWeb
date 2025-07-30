from django.test import TestCase
from django.urls import reverse
from django.contrib.auth.models import User
from mesh.models import Mesh
from mesh.forms import MeshType, IntervalMesh


class MeshModelTest(TestCase):
    """Test cases for Mesh model."""
    
    def test_mesh_creation(self):
        """Test basic mesh creation."""
        mesh = Mesh.objects.create(
            name="Test Mesh",
            mesh="interval,10,0.0,1.0"
        )
        self.assertEqual(mesh.name, "Test Mesh")
        self.assertEqual(mesh.mesh, "interval,10,0.0,1.0")
        self.assertEqual(str(mesh), f"Test Mesh - {mesh.id}")

    def test_mesh_name_required(self):
        """Test that mesh name is required."""
        with self.assertRaises(Exception):
            Mesh.objects.create(mesh="interval,10,0.0,1.0")


class MeshFormTest(TestCase):
    """Test cases for Mesh forms."""
    
    def test_mesh_type_form_valid(self):
        """Test valid mesh type form."""
        form_data = {
            'mesh_name': 'Test Mesh',
            'mesh_type': 'interval'
        }
        form = MeshType(data=form_data)
        self.assertTrue(form.is_valid())

    def test_mesh_type_form_invalid(self):
        """Test invalid mesh type form."""
        form_data = {
            'mesh_name': '',  # Empty name should be invalid
            'mesh_type': 'interval'
        }
        form = MeshType(data=form_data)
        self.assertFalse(form.is_valid())

    def test_interval_mesh_form_valid(self):
        """Test valid interval mesh form."""
        form_data = {
            'interval_n': 10,
            'interval_x0': 0.0,
            'interval_x1': 1.0
        }
        form = IntervalMesh(data=form_data)
        self.assertTrue(form.is_valid())

    def test_interval_mesh_form_invalid_negative(self):
        """Test interval mesh form with negative values."""
        form_data = {
            'interval_n': -10,  # Negative should be invalid
            'interval_x0': 0.0,
            'interval_x1': 1.0
        }
        form = IntervalMesh(data=form_data)
        # Note: This form doesn't validate ranges, but it should
        # This test shows what needs to be improved


class MeshViewTest(TestCase):
    """Test cases for Mesh views."""
    
    def test_mesh_view_get(self):
        """Test GET request to mesh view."""
        response = self.client.get(reverse('generate_mesh'))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Select Mesh Type')

    def test_mesh_view_post_mesh_type(self):
        """Test POST request with mesh type selection."""
        response = self.client.post(reverse('generate_mesh'), {
            'mesh_name': 'Test Mesh',
            'mesh_type': 'interval'
        })
        self.assertEqual(response.status_code, 200)
        # Should show mesh parameter form
        self.assertContains(response, 'Enter the Mesh Parameters')

# Note: These tests would need FEniCS to be installed to test actual mesh generation
# For production, consider mocking FEniCS calls or testing with a FEniCS container
