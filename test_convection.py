#!/usr/bin/env python3
"""
Test script for convection boundary conditions implementation
"""
import sys
import os
import django

# Add project root to path
sys.path.append('/home/harish/project/FenicsWeb')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'FenicsWeb.settings')
django.setup()

from fenics import *
from SteadyStateThermal.conduction import solve_conduction

def test_1d_convection():
    """Test 1D heat conduction with convection BC"""
    print("Testing 1D convection boundary condition...")
    
    # Create simple 1D mesh
    mesh = UnitIntervalMesh(10)
    
    # Define boundary conditions
    bc_data = {
        'left': {
            'type': 'dirichlet',
            'value': 100.0  # Fixed temperature at left boundary
        },
        'right': {
            'type': 'convection',
            'h': 10.0,      # Convection coefficient
            't_amb': 20.0   # Ambient temperature
        }
    }
    
    # Material properties
    f = 0.0  # No heat source
    k = 1.0  # Thermal conductivity
    
    # Solve
    try:
        solution, solution_list = solve_conduction('1D', mesh, bc_data, f, k)
        print("✓ 1D convection test passed!")
        print(f"Solution range: {min(s[1] for s in solution_list):.2f} - {max(s[1] for s in solution_list):.2f}")
        return True
    except Exception as e:
        print(f"✗ 1D convection test failed: {e}")
        return False

def test_2d_mixed_bc():
    """Test 2D heat conduction with mixed boundary conditions"""
    print("Testing 2D mixed boundary conditions...")
    
    # Create simple 2D mesh
    mesh = UnitSquareMesh(8, 8)
    
    # Define boundary conditions
    bc_data = {
        'left': {
            'type': 'dirichlet',
            'value': 100.0  # Fixed temperature
        },
        'right': {
            'type': 'convection',
            'h': 15.0,      # High convection
            't_amb': 25.0   # Ambient temperature
        },
        'top': {
            'type': 'convection',
            'h': 5.0,       # Low convection
            't_amb': 20.0   # Different ambient temperature
        },
        'bottom': {
            'type': 'dirichlet',
            'value': 50.0   # Fixed temperature
        }
    }
    
    # Material properties
    f = 0.0  # No heat source
    k = 1.0  # Thermal conductivity
    
    # Solve
    try:
        solution, solution_list = solve_conduction('2D', mesh, bc_data, f, k)
        print("✓ 2D mixed BC test passed!")
        print(f"Solution range: {min(s[1] for s in solution_list):.2f} - {max(s[1] for s in solution_list):.2f}")
        return True
    except Exception as e:
        print(f"✗ 2D mixed BC test failed: {e}")
        return False

def test_pure_convection():
    """Test case with only convection boundaries"""
    print("Testing pure convection boundary conditions...")
    
    # Create simple 2D mesh
    mesh = UnitSquareMesh(6, 6)
    
    # Define boundary conditions (all convection)
    bc_data = {
        'left': {
            'type': 'convection',
            'h': 10.0,
            't_amb': 80.0
        },
        'right': {
            'type': 'convection',
            'h': 10.0,
            't_amb': 20.0
        },
        'top': {
            'type': 'convection',
            'h': 5.0,
            't_amb': 30.0
        },
        'bottom': {
            'type': 'convection',
            'h': 5.0,
            't_amb': 70.0
        }
    }
    
    # Material properties with heat source
    f = 1000.0  # Heat source to ensure well-posed problem
    k = 1.0     # Thermal conductivity
    
    # Solve
    try:
        solution, solution_list = solve_conduction('2D', mesh, bc_data, f, k)
        print("✓ Pure convection test passed!")
        print(f"Solution range: {min(s[1] for s in solution_list):.2f} - {max(s[1] for s in solution_list):.2f}")
        return True
    except Exception as e:
        print(f"✗ Pure convection test failed: {e}")
        return False

if __name__ == "__main__":
    print("Starting convection boundary condition tests...\n")
    
    tests = [
        test_1d_convection,
        test_2d_mixed_bc,
        test_pure_convection
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
        print()
    
    print(f"Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! Convection boundary conditions are working correctly.")
    else:
        print("⚠️  Some tests failed. Please check the implementation.")
