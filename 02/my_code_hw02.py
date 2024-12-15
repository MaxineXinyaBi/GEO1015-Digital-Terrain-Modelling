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
        triangle_idx = dt.locate([x,y])
        # step2: the coordinates of triangle vertices
        p1 = dt.get_point(triangle_idx[0])
        p2 = dt.get_point(triangle_idx[1])
        p3 = dt.get_point(triangle_idx[2])
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


def interpolate_laplace(dt, x, y):
    """
    !!! TO BE COMPLETED !!!
    !!! You are free to subdivide the functionality of this function into several functions !!!

    Function that interpolates at location (x,y) in a DT with the Laplace interpolation.

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
    #return dt.interpolate({"method": "Laplace"}, [[x, y]])
    # -- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    try:
        # step1: find the triangle that has the point
        triangle_idx = dt.locate([x,y])
        v1 = dt.get_point(triangle_idx[0])
        v2 = dt.get_point(triangle_idx[1])
        v3 = dt.get_point(triangle_idx[2])
        x1, y1, z1 = v1
        x2, y2, z2 = v2
        x3, y3, z3 = v3
        vertices = [v1, v2, v3]
        edges = [(v1, v2), (v2, v3), (v3, v1)]
        # step2: the unknown point is on the vertice, return the z value of the vertice
        allowed_error_laplace = 1e-10
        if abs(x - x1) <= allowed_error_laplace and abs(y - y1) <= allowed_error_laplace:
            return z1
        if abs(x - x2) <= allowed_error_laplace and abs(y - y2) <= allowed_error_laplace:
            return z2
        if abs(x - x3) <= allowed_error_laplace and abs(y - y3) <= allowed_error_laplace:
            return z3
        # step3: the edge in voronoi -- the height from point to another edge，
        # step4: the distance between point to the centre of the cell -- distance between unknown point to opposite vertice
        total_area = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2
        weights = []
        z_value_vertice = []

        for i, edge in enumerate(edges):
            edge_length = math.sqrt((edge[0][0] - edge[1][0]) ** 2 + (edge[0][1] - edge[1][1]) ** 2)
            height = 2 * total_area / edge_length
            opposite_vertice = vertices[(i + 2) % 3]
            dist_to_vertice = math.sqrt((x - opposite_vertice[0])**2 + (y - opposite_vertice[1])**2)
            weight = dist_to_vertice / height
            weights.append(weight)
            z_value_vertice.append(opposite_vertice[2])

        # step5: calculate unknown z
        unknown_z = sum(w * z for w, z in zip(weights, z_value_vertice)) / sum(weights)
        return unknown_z
    except Exception as e:
        return np.nan



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
    pass
    # # -- generate 100 points randomly in the plane
    # rng = np.random.default_rng(seed=42)
    # pts = rng.random((100, 3))
    # pts = pts * 100
    # dt = startinpy.DT()
    # dt.insert(pts, insertionstrategy="AsIs")
    # # -- showcase for tqdm package to see progress of a loop
    # for i in tqdm.tqdm(range(500)):
    #     time.sleep(0.01)
    # return dt

def filter_outlier(pts, k = 1.5):
    """the function is to filter out outlier points using IQR method"""
    z_values = pts[:, 2]
    q1 = np.percentile(z_values, 25)
    q3 = np.percentile(z_values, 75)
    iqr = q3 - q1
    interquartile_range =  iqr * k

    lower_fence = q1 - interquartile_range
    upper_fence = q3 + interquartile_range

    mask = (z_values >= lower_fence) & (z_values <= upper_fence)
    filtered_pts = pts[mask]
    outliers = pts[~mask]
    return filtered_pts, outliers

def find_lowest_pt_in_cell(filter_pts, resolution):
    """create initial TIN using the lowest pt in each cell"""
    x_min, x_max = np.min(filter_pts[:, 0]), np.max(filter_pts[:, 0])
    y_min, y_max = np.min(filter_pts[:, 1]), np.max(filter_pts[:, 1])

    # calculate the number of cells
    x_cell = np.ceil((x_max - x_min) / resolution).astype(int)
    y_cell = np.ceil((y_max - y_min) / resolution).astype(int)

    # loop every cell
    initial_tin_pts = []
    for i in range(x_cell):
        for j in range(y_cell):
            x_left = x_min + i * resolution
            x_right = x_left + (i + 1) * resolution
            y_bottom = y_min + j * resolution
            y_top = y_min + (j + 1) * resolution

            # create the mask to filter out the points that are in the current cell
            pts_in_cell_mask = (filter_pts[:,0] >= x_left) & (filter_pts[:,0] < x_right) & (filter_pts[:,1] >= y_bottom) & (filter_pts[:,1] < y_top)
            pts_in_cell = filter_pts[pts_in_cell_mask]
            # if there are points in the cell find the point with lowest z value and add it to the initial tin vertices list
            if len(pts_in_cell) > 0:
                lowest_pt_index = np.argmin(pts_in_cell[:, 2])
                initial_tin_pts.append(pts_in_cell[lowest_pt_index])

    return np.array(initial_tin_pts)

def create_initial_tin(qualified_pts):
    """create initial tin with the lowest pt in each cell"""
    initial_dt = startinpy.DT()
    initial_dt.insert(qualified_pts)
    return initial_dt

def remaining_pts(filtered_pts, initial_tin_pts):
    """get the remaining pt"""
    mask = ~np.isin(filtered_pts, initial_tin_pts).all(axis=1)
    remaining_pts = filtered_pts[mask]
    return remaining_pts

def densification_tin(initial_dt,remaining_pts,max_dist,max_angle):
    pass
