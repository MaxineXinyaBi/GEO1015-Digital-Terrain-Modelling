import numpy as np
import startinpy
import my_code_hw02
import math


def compare_with_startinpy(dt, test_points, test_name):
    """Compare our implementation with startinpy's built-in interpolation"""
    print(f"\n=== {test_name} ===")

    for x, y in test_points:
        # Our implementations
        our_linear = my_code_hw02.interpolate_linear(dt, x, y)
        our_laplace = my_code_hw02.interpolate_laplace(dt, x, y)

        # Startinpy implementations (for comparison only)
        startinpy_linear = dt.interpolate({"method": "TIN"}, [[x, y]])[0]
        startinpy_laplace = dt.interpolate({"method": "Laplace"}, [[x, y]])[0]

        print(f"\nTest point ({x}, {y}):")
        print(f"Linear - Our implementation: {our_linear:.6f}")
        print(f"Linear - Startinpy: {startinpy_linear:.6f}")
        print(f"Difference: {abs(our_linear - startinpy_linear):.6f}")

        print(f"Laplace - Our implementation: {our_laplace:.6f}")
        print(f"Laplace - Startinpy: {startinpy_laplace:.6f}")
        print(f"Difference: {abs(our_laplace - startinpy_laplace):.6f}")


def run_tests():
    # Create a simple triangle mesh
    pts = np.array([
        [0, 0, 0],  # Point 1
        [3, 0, 3],  # Point 2
        [0, 3, 6],  # Point 3
        [3, 3, 9]  # Point 4 (creates two triangles)
    ])

    dt = startinpy.DT()
    dt.insert(pts)

    # Test cases
    vertex_test = [(0, 0)]  # Point on vertex
    edge_test = [(1.5, 0)]  # Point on edge
    inside_test = [(1, 1)]  # Point inside triangle
    outside_test = [(5, 5)]  # Point outside convex hull

    # Run comparisons
    compare_with_startinpy(dt, vertex_test, "Point on vertex")
    compare_with_startinpy(dt, edge_test, "Point on edge")
    compare_with_startinpy(dt, inside_test, "Point inside triangle")
    compare_with_startinpy(dt, outside_test, "Point outside convex hull")


if __name__ == "__main__":
    print("Starting comparison tests...")
    run_tests()