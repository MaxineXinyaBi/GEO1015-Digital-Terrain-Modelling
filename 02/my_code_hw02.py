# -- my_code_hw02.py
# -- geo1015.2024.hw02
# -- [YOUR NAME]
# -- [YOUR STUDENT NUMBER]


import time
import math
import numpy as np
import startinpy
import tqdm


def interpolate_linear(dt, x, y):
    """
    !!! TO BE COMPLETED !!!
    !!! You are free to subdivide the functionality of this function into several functions !!!

    Function that interpolates at location (x,y) in a DT with the linear in TIN interpolation.

    Inputs:
      dt: the startinpy DT
      x:  coordinate x to interpolate
      y:  coordinate y to interpolate

    Output:
      - the estimated value for z
      - np.nan if outside the convex hull (impossible to interpolate)
        (NaN: https://numpy.org/doc/stable/reference/constants.html#numpy.nan)
    """
    # -- !!! dt.interpolate() is illegal to use for this assignment
    # -- !!! you need to write your own code, and you can use the functions in startinpy
    # return dt.interpolate({"method": "TIN"}, [[x, y]])
    # -- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    try:
        # step1: find the triangle that has the point
        containing_triangle = dt.locate([x,y])
        if containing_triangle is not None and dt.is_finite(containing_triangle):
            # step2: the coordinates of triangle vertices
            p1 = dt.get_point(containing_triangle[0])
            p2 = dt.get_point(containing_triangle[1])
            p3 = dt.get_point(containing_triangle[2])
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            x3, y3, z3 = p3

        # step3: if the point is on the vertex, return the z value of the vertex
        allowed_error_linear = 1e-10
        if abs(x - x1) <= allowed_error_linear and abs(y - y1) <= allowed_error_linear:
            return z1
        if abs(x - x2) <= allowed_error_linear and abs(y - y2) <= allowed_error_linear:
            return z2
        if abs(x - x3) <= allowed_error_linear and abs(y - y3) <= allowed_error_linear:
            return z3

        # step4: if the point is on the edge/inside of the triangle
        total_area = abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)) / 2
        area1 = abs((x2 - x) * (y3 - y) - (x3 - x) * (y2 - y)) / 2
        area2 = abs((x3 - x) * (y1 - y) - (x1 - x) * (y3 - y)) / 2
        area3 = abs((x1 - x) * (y2 - y) - (x2 - x) * (y1 - y)) / 2

        w1 = area1 / total_area
        w2 = area2 / total_area
        w3 = area3 / total_area

        unknown_z = w1 * z1 + w2 * z2 + w3 * z3

        return unknown_z
    except Exception as e:
        print(e)
        return np.nan


def find_natural_neighbors(dt, x, y):
    temp_point = [x, y, 0.0]
    result = dt.insert_one_pt(temp_point)

    temp_point_idx = result[0]
    is_new_vertex = result[1]

    if not is_new_vertex:
        z_value = dt.get_point(temp_point_idx)[2]
        dt.remove(temp_point_idx)
        return z_value

    try:
        neighbors = dt.adjacent_vertices_to_vertex(temp_point_idx)
    except:
        dt.remove(temp_point_idx)
        return None

    if len(neighbors ) == 0:
        dt.remove(temp_point_idx)
        return None

    dt.remove(temp_point_idx)
    return neighbors

def get_circumcircle(p1, p2, p3):
    """calculate circumcircle from triangle"""
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]
    x3, y3 = p3[0], p3[1]

    d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < 1e-10:
        return None, None

    ux = ((x1 * x1 + y1 * y1) * (y2 - y3) +
          (x2 * x2 + y2 * y2) * (y3 - y1) +
          (x3 * x3 + y3 * y3) * (y1 - y2)) / d
    uy = ((x1 * x1 + y1 * y1) * (x3 - x2) +
          (x2 * x2 + y2 * y2) * (x1 - x3) +
          (x3 * x3 + y3 * y3) * (x2 - x1)) / d

    center = (ux, uy)
    radius = ((center[0] - x1) ** 2 + (center[1] - y1) ** 2) ** 0.5

    return center, radius


def get_circumcenters(dt, point, neighbors):
    centers = []
    n = len(neighbors)

    for i in range(n):
        p1 = dt.get_point(neighbors[i])
        p2 = dt.get_point(neighbors[(i + 1) % n])

        center, _ = get_circumcircle(p1, p2, point)
        if center is not None:
            centers.append(center)

    return centers


def calculate_weights(centers, neighbors, dt, x, y):
    """calculate weights based on centers and neighbors"""
    weights = []
    for i in range(len(neighbors)):
        center1 = centers[i]
        center2 = centers[i - 1]
        voronoi_length = np.sqrt(
            (center1[0] - center2[0]) ** 2 +(center1[1] - center2[1]) ** 2)

        neighbor_point = dt.get_point(neighbors[i])
        distance = np.sqrt((x - neighbor_point[0]) ** 2 + (y - neighbor_point[1]) ** 2)

        if distance > 1e-10:
            weight = voronoi_length / distance
            weights.append({
                'weight': weight,
                'z': neighbor_point[2]
            })

    return weights


def calculate_interpolation(weights):
    total_weight = sum(w['weight'] for w in weights)
    if total_weight <= 0:
        return None

    unknown_z = sum(w['weight'] * w['z'] for w in weights) / total_weight
    return unknown_z


def is_on_voronoi_edge(center1, center2, point, tolerance=1e-10):
    """point on voronoi edge"""
    x1, y1 = center1
    x2, y2 = center2
    x, y = point

    if abs(x1 - x2) < tolerance and abs(y1 - y2) < tolerance:
        return False

    dx = x2 - x1
    dy = y2 - y1
    t = ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy)

    if t < 0 or t > 1:
        return False

    px = x1 + t * dx
    py = y1 + t * dy

    return abs(x - px) < tolerance and abs(y - py) < tolerance


def interpolate_laplace(dt, x, y):
    """Laplace interpolation at point (x,y)"""
    # Check if point is inside convex hull
    try:
        containing_triangle = dt.locate([x, y])
        if containing_triangle is None or not dt.is_finite(containing_triangle):
            return np.nan
    except:
        return np.nan

    # check if the point overlay with known points
    result = find_natural_neighbors(dt, x, y)
    if isinstance(result, (float, np.float64)):
        return result
    if result is None:
        return np.nan
    neighbors = result

    # circumcenter
    centers = get_circumcenters(dt, [x, y, 0.0], neighbors)
    if len(centers) < 3:
        return np.nan

    # point on voronoi edge
    for i in range(len(centers)):
        center1 = centers[i]
        center2 = centers[(i + 1) % len(centers)]

        if is_on_voronoi_edge(center1, center2, [x, y]):
            p1 = dt.get_point(neighbors[i])
            p2 = dt.get_point(neighbors[(i + 1) % len(neighbors)])
            d1 = np.sqrt((x - p1[0]) ** 2 + (y - p1[1]) ** 2)
            d2 = np.sqrt((x - p2[0]) ** 2 + (y - p2[1]) ** 2)
            total_d = d1 + d2
            w1 = d2 / total_d
            w2 = d1 / total_d
            return w1 * p1[2] + w2 * p2[2]

    # weights and final results
    weights = calculate_weights(centers, neighbors, dt, x, y)

    unknown_z = calculate_interpolation(weights)
    if unknown_z is None:
        return np.nan
    else:
        return unknown_z


def gftin(pts, resolution, max_dist, max_angle):
    """
    !!! TO BE COMPLETED !!!
    !!! the code written below is just dummy code, replace it !!!
    !!! You are free to subdivide the functionality of this function into several functions !!!

    Function that performs ground filtering using TIN refinement and returns the DT.

    Inputs:
      pts:        the Nx3 numpy array of the PC
      resolution: resolution (cellsize) for the initial grid that is computed as part of the ground filtering algorithm,
      max_dist:   distance threshold used in the ground filtering algorithm,
      max_angle:  angle threshold used in the ground filtering algorithm in degrees,

    Output:
      the startinpy DT of the ground
    """
    # # step1: filter outlier
    # filtered_points = filter_outlier(pts, k = 1.5)

    # step2 : construct initial tin network
    initial_tin_vertices = find_lowest_pt_in_cell_with_bbox(pts, resolution)
    initial_tin = create_initial_tin(initial_tin_vertices)


    # step3: get remaining points and construct final tin network
    remaining_points = remaining_pts(pts, initial_tin_vertices)
    final_tin = densification_tin(initial_tin, remaining_points, max_dist, max_angle)

    return final_tin


# def filter_outlier(pts, k = 1.5):
#     """the function is to filter out outlier points using IQR method"""
#     z_values = pts[:, 2]
#     q1 = np.percentile(z_values, 25)
#     q3 = np.percentile(z_values, 75)
#     iqr = q3 - q1
#     interquartile_range =  iqr * k
#
#     lower_fence = q1 - interquartile_range
#     upper_fence = q3 + interquartile_range
#
#     mask = (z_values >= lower_fence) & (z_values <= upper_fence)
#     filtered_pts = pts[mask]
#     return filtered_pts


def find_lowest_pt_in_cell_with_bbox(pts, resolution):
    """Create initial TIN using the lowest point in each cell and add bounding box corner points."""
    # Convert to numpy array (if not already)
    pts = np.array(pts)

    # Get bounding box
    x_min, y_min = np.min(pts[:, 0]), np.min(pts[:, 1])
    x_max, y_max = np.max(pts[:, 0]), np.max(pts[:, 1])
    bbox = [x_min, y_min, x_max, y_max]
    print(f"Bounding box: {bbox}")

    # Calculate the number of cells
    x_cell = int(np.ceil((x_max - x_min) / resolution))
    y_cell = int(np.ceil((y_max - y_min) / resolution))

    # Loop through every cell to find the lowest point
    initial_tin_pts = []
    for i in range(x_cell):
        for j in range(y_cell):
            x_left = x_min + i * resolution
            x_right = x_min + (i + 1) * resolution
            y_bottom = y_min + j * resolution
            y_top = y_min + (j + 1) * resolution

            # Create the mask to filter points in the current cell
            pts_in_cell_mask = (pts[:, 0] >= x_left) & (pts[:, 0] < x_right) & (pts[:, 1] >= y_bottom) & (pts[:, 1] < y_top)
            pts_in_cell = pts[pts_in_cell_mask]

            # If there are points in the cell, find the point with the lowest z value
            if len(pts_in_cell) > 0:
                lowest_pt_index = np.argmin(pts_in_cell[:, 2])
                initial_tin_pts.append(pts_in_cell[lowest_pt_index])

    # Convert to numpy array
    initial_tin_pts = np.array(initial_tin_pts)

    # Define bounding box corner points (2D)
    buffer = 0.5
    bbox_corners = np.array([
        [x_min- buffer, y_min - buffer],
        [x_max + buffer, y_min - buffer],
        [x_max + buffer, y_max + buffer],
        [x_min - buffer, y_max + buffer],
    ])

    # Find the closest point in the point cloud to each corner
    for corner in bbox_corners:
        distances = np.linalg.norm(pts[:, :2] - corner, axis=1)
        closest_point_index = np.argmin(distances)
        closest_point = pts[closest_point_index]  # Get the 3D coordinate [x, y, z]
        print(f"Corner: {corner}, Closest Point: {closest_point}, Distance: {distances[closest_point_index]}")
        initial_tin_pts = np.vstack([initial_tin_pts, closest_point])  # Add directly to the TIN points

    return initial_tin_pts


def create_initial_tin(qualified_pts):
    """create initial tin with the lowest pt in each cell"""
    initial_dt = startinpy.DT()
    initial_dt.insert(qualified_pts)
    return initial_dt

def remaining_pts(pts, initial_tin_pts):
    """Get the remaining points"""
    mask = ~np.any(np.all(pts[:, None, :] == initial_tin_pts[None, :, :], axis=2), axis=1)
    return pts[mask]

def densification_tin(initial_dt, remaining_pts, max_dist, max_angle):
    """the densification of initial tin network"""

    print(f"Starting points to process: {len(remaining_pts)}")
    remaining_pts_list = remaining_pts.tolist()
    outside_convex_hull_points = []  # 用于记录凸包外的点
    iteration = 0

    while remaining_pts_list:
        iteration += 1
        points_added = False
        points_checked = 0
        points_added_this_round = 0

        for pt in remaining_pts_list[:]:
            points_checked += 1
            triangle_idx = initial_dt.locate([pt[0], pt[1]])

            if triangle_idx is None:
                # 如果 locate 返回 None，说明点在凸包之外
                outside_convex_hull_points.append(pt)
                print(f"Point outside convex hull: x={pt[0]}, y={pt[1]}, z={pt[2]}")
                raise ValueError("Point outside convex hull")
                continue

            p1 = initial_dt.get_point(triangle_idx[0])
            p2 = initial_dt.get_point(triangle_idx[1])
            p3 = initial_dt.get_point(triangle_idx[2])

            distance = pt_distance_to_plane(pt, p1, p2, p3)
            angle = calculate_max_angle(pt, p1, p2, p3)

            if distance < max_dist and angle < max_angle:
                initial_dt.insert([pt])
                remaining_pts_list.remove(pt)
                points_added = True
                points_added_this_round += 1

        print(f"Iteration {iteration}:")
        print(f"  Points checked: {points_checked}")
        print(f"  Points added: {points_added_this_round}")
        print(f"  Remaining points: {len(remaining_pts_list)}")

        if not points_added:
            print("No points added in this iteration, breaking")
            break

    # 打印凸包外点的数量
    print(f"Total points outside convex hull: {len(outside_convex_hull_points)}")

    # 如果有点在凸包外，打印这些点的坐标
    if outside_convex_hull_points:
        print("Points outside convex hull:")
        for pt in outside_convex_hull_points:
            print(f"  x: {pt[0]:.3f}, y: {pt[1]:.3f}, z: {pt[2]:.3f}")

    return initial_dt




def calculate_normal_vector(p1, p2, p3):
    """calculate the normal vector of a plane"""
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    normal = normal / np.linalg.norm(normal)  # normalize normal vector of the plane
    return normal


def pt_distance_to_plane(pt,p1,p2,p3):
    """calculate the distance between point and plane"""
    normal_vec = calculate_normal_vector(p1,p2,p3)
    ptp1 = pt - p1
    distance = abs(np.dot(ptp1,normal_vec))
    return distance


def calculate_max_angle(pt, p1, p2, p3):
    """calculate the max angle of pt to p1 p2 p3"""
    normal_vec = calculate_normal_vector(p1,p2,p3)
    dist = pt_distance_to_plane(pt,p1,p2,p3)
    proj_pt = pt - normal_vec * dist

    angles = []
    for vertex in [p1, p2, p3]:
        proj_pt_to_vertex = np.linalg.norm(proj_pt - vertex)
        degree = np.degrees(np.arctan2( dist,proj_pt_to_vertex))
        angles.append(float(degree))

    return max(angles)




