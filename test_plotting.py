#!/usr/bin/env python3
"""
Test script for the plotting function
"""
import sys
import os
import django

# Add project root to path
sys.path.append('/home/harish/project/FenicsWeb')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FenicsWeb.settings')
django.setup()

from fenics import *
from SteadyStateThermal.conduction import solve_conduction, generate_solution_plot

def test_plotting():
    """Test the plotting function"""
    print("Testing plotting function...")
    
    # Create simple 1D mesh
    mesh = UnitIntervalMesh(10)
    
    # Define boundary conditions
    bc_data = {
        'left': {
            'type': 'dirichlet',
            'value': 100.0
        },
        'right': {
            'type': 'convection',
            'h': 10.0,
            't_amb': 20.0
        }
    }
    
    # Solve
    try:
        solution, solution_list = solve_conduction('1D', mesh, bc_data, 0.0, 1.0)
        plot_filename = generate_solution_plot(solution)
        print(f"✓ 1D plotting test passed! Plot saved as: {plot_filename}")
        
        # Test 2D plotting
        mesh_2d = UnitSquareMesh(8, 8)
        bc_data_2d = {
            'left': {'type': 'dirichlet', 'value': 100.0},
            'right': {'type': 'convection', 'h': 15.0, 't_amb': 25.0},
            'top': {'type': 'convection', 'h': 5.0, 't_amb': 20.0},
            'bottom': {'type': 'dirichlet', 'value': 50.0}
        }
        
        solution_2d, _ = solve_conduction('2D', mesh_2d, bc_data_2d, 0.0, 1.0)
        plot_filename_2d = generate_solution_plot(solution_2d)
        print(f"✓ 2D plotting test passed! Plot saved as: {plot_filename_2d}")
        
        return True
    except Exception as e:
        print(f"✗ Plotting test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    test_plotting()
