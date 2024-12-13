# -- geo1015.2024.hw01
# -- [YOUR NAME] Xinya Bi
# -- [YOUR STUDENT NUMBER] 6195350

import argparse
import json
import math
import sys
import rasterio
import startinpy
import numpy as np


def sample_raster(dataset, thining):
    """The function use the input geotiff file, check no data value, perform random sampling based on required scale"""
    try:
        with rasterio.open(dataset) as d:
            # print input file info
            print("Dataset Info:")
            print(f"Driver: {d.driver}")
            print(f"File: {d.name}")
            print(f"Size: {d.width} x {d.height} pixels")
            print(f"Bounds: {d.bounds}")
            print(f"CRS: {d.crs}")
            print(f"Resolution: {d.res}")
            print(f"No data value: {d.nodata}")

            if d.count != 1:
                raise ValueError("Input must be a single-band GeoTIFF")

            # read data and get valid data
            raster_data = d.read(1, masked=True)
            valid_values = raster_data.compressed()

            # calculate random sample num
            sample_num = int(len(valid_values) * thining)

            # random sampling based on sample_num
            sample_pixels = np.random.choice(len(valid_values),sample_num,replace=False)

            # extract pos and value of sample points
            valid_positions = np.where(~raster_data.mask)
            sample_values = valid_values[sample_pixels]
            sample_row = valid_positions[0][sample_pixels]
            sample_col = valid_positions[1][sample_pixels]

            # transform to actual coords
            sample_pixels_coords = []
            for i, (row, column) in enumerate(zip(sample_row, sample_col)):
                x, y = d.transform * (column, row)
                sample_pixels_coords.append((x, y, sample_values[i]))

            return sample_pixels_coords

    except rasterio.errors.RasterioIOError as e:
        print(f"Error opening file: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Invalid data: {e}")
        sys.exit(1)

def create_tin(points):
    """create tin network with all random sampled points"""
    dt = startinpy.DT()
    dt.insert(points)

    # get triangles, vertices(except infinite point) convex hull
    triangles = []
    for triangle in dt.triangles:
        if dt.is_finite(triangle):
            triangles.append(triangle)
    vertices = dt.points[1:]
    hull = dt.convex_hull()

    return dt, triangles, vertices, hull


def get_intersection(dt, triangle, height, tolerance=1e-10):
    """return the intersection point of a triangle and contour line/or no intersection(none)"""
    points = dt.points[1:]
    v1 = points[triangle[0]-1]
    v2 = points[triangle[1]-1]
    v3 = points[triangle[2]-1]

    # case1: check if height is in the range of triangle, if not return none
    z_min = min(v1[2],v2[2],v3[2]) - tolerance
    z_max = max(v1[2],v2[2],v3[2]) + tolerance
    if height < z_min or height > z_max:
        return None

    intersections = []
    vertices =[v1,v2,v3]
    edges = [(v1,v2),(v2,v3),(v3,v1)]

    # deal with flat triangles
    # check flat triangles (if they are on contour lines)
    z_diffs = [abs(v1[2] - v2[2]),abs(v2[2] - v3[2]),abs(v3[2] - v1[2])]
    num_flat_edges = 0
    for diff in z_diffs:
        if diff <= tolerance and abs(v1[2] - height) < tolerance and abs(v2[2] - height) < tolerance and abs(v3[2] - height) < tolerance:
            num_flat_edges += 1
    # case2: complete flat triangles, add no vertices
    if num_flat_edges == 3:
        return None
    # case3: triangle has 1 or 2 flat edges, and they are on contour lines
    elif num_flat_edges == 1:
        # add intersection based on normal vectors
        for edge in edges:
            v1, v2 = edge
            if abs(v1[2] - v2[2]) <= tolerance:
                dx = v2[0] - v1[0]
                dy = v2[1] - v1[1]
                normal_vector_x = -dy
                normal_vector_y = dx

                # if normal vector's y > 0 or (y = 0 and x > 0) then add
                if (normal_vector_y > tolerance) or (abs(normal_vector_y) < tolerance and normal_vector_x > tolerance):
                    intersections.append((v1[0], v1[1], height))
                    intersections.append((v2[0], v2[1], height))

    else:
        # if the triangle doesn't have flat edges
        # case4: contour cross the vertice of triangle and an edge of triangle
        # case5: contour cross the 2 edges of triangle
        for vertice in vertices:
            if abs(vertice[2] - height) <= tolerance:
                intersections.append((vertice[0], vertice[1], height))

        # check for height cross triangle edge
        for edge in edges:
            vertex1, vertex2 = edge

            if (vertex1[2] - height) * (vertex2[2] - height) < 0:
                if vertex1[2] > vertex2[2]:
                    high_end = vertex1
                    low_end = vertex2
                else:
                    high_end = vertex2
                    low_end = vertex1

                dz = high_end[2] - low_end[2]
                if abs(dz) > tolerance:
                    linescale = (height - low_end[2]) / dz
                    intersect_x = low_end[0] + linescale * (high_end[0] - low_end[0])
                    intersect_y = low_end[1] + linescale * (high_end[1] - low_end[1])
                    intersections.append((intersect_x, intersect_y, height))

    return intersections

def get_all_intersections(dt, triangles, heights, tolerance=1e-10):
    intersections_by_height = {height:[] for height in heights}
    points = dt.points[1:]
    for triangle in triangles:
        v1 = points[triangle[0]-1]
        v2 = points[triangle[1]-1]
        v3 = points[triangle[2]-1]

        z_min = min(v1[2], v2[2], v3[2]) - tolerance
        z_max = max(v1[2], v2[2], v3[2]) + tolerance

        for height in heights:
            if height < z_min or height > z_max:
                continue

            intersections = get_intersection(dt, triangle, height, tolerance=tolerance)
            if intersections:
                formatted_intersections = [(round(float(x), 6), round(float(y), 6), round(float(z), 6))for x, y,z in intersections]
                intersections_by_height[height].append(formatted_intersections)
    return intersections_by_height


def structure_contours(contour_dict, dt):
    """the left side is always higher than right side"""
    result = {}
    points = dt.points[1:]

    for height, segments in contour_dict.items():
        # filter out single point
        valid_segments = [seg for seg in segments if len(seg) > 1]
        all_ordered_segments = []

        for segment in valid_segments:
            # get the mid point of segment
            mid_x = (segment[0][0] + segment[1][0]) / 2
            mid_y = (segment[0][1] + segment[1][1]) / 2

            containing_triangle = None
            for triangle in dt.triangles:
                if dt.is_finite(triangle):
                    v1 = points[triangle[0] - 1]
                    v2 = points[triangle[1] - 1]
                    v3 = points[triangle[2] - 1]

                    if is_point_in_triangle((mid_x, mid_y),
                                            (v1[0], v1[1]),
                                            (v2[0], v2[1]),
                                            (v3[0], v3[1])):
                        containing_triangle = (v1, v2, v3)
                        break

            if containing_triangle:
                dx = segment[1][0] - segment[0][0]
                dy = segment[1][1] - segment[0][1]

                # calculate the point that's on the left side of contour
                epsilon = 1e-6
                left_x = mid_x - dy * epsilon
                left_y = mid_y + dx * epsilon
                left_z = point_in_triangle_height((left_x, left_y), containing_triangle)

                # if left side is lower, flip the edge
                if left_z < height:
                    segment = segment[::-1]

            all_ordered_segments.append(segment)

        result[height] = all_ordered_segments

    return result

def sign(p1, p2, p3):
    """ if > 0, p1 is on the left side of p2p3"""
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

def is_point_in_triangle(p, v1, v2, v3):
    d1 = sign(p, v1, v2)
    d2 = sign(p, v2, v3)
    d3 = sign(p, v3, v1)

    has_negative = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_positive = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not (has_negative and has_positive)

def area(x1, y1, x2, y2, x3, y3):
    """triangle area"""
    return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)


def point_in_triangle_height(p, triangle):
    """Barycentric coordinate"""
    v1, v2, v3 = triangle

    A = area(v1[0], v1[1], v2[0], v2[1], v3[0], v3[1])

    A1 = area(p[0], p[1], v2[0], v2[1], v3[0], v3[1])
    A2 = area(v1[0], v1[1], p[0], p[1], v3[0], v3[1])
    A3 = area(v1[0], v1[1], v2[0], v2[1], p[0], p[1])

    w1 = A1 / A
    w2 = A2 / A
    w3 = A3 / A

    return w1 * v1[2] + w2 * v2[2] + w3 * v3[2]


def calculate_azimuth(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    dx = x2 - x1
    dy = y2 - y1
    azimuth = 180 / math.pi * math.atan2(dx, dy)
    return azimuth % 360

def calculate_lightness(azimuth):
    l = abs(((azimuth - 225) % 360) - 180)
    lightness = l / 180 * 100
    return lightness

def create_json(ordered_contours):
    mygeojson = {}
    mygeojson['type'] = 'FeatureCollection'
    mygeojson["features"] = []
    for height, edges in ordered_contours.items():
        for edge in edges:
            # only keep complete lines
            if len(edge) != 1:
                f = {}
                f["type"] = "Feature"
                start_point = [edge[0][0], edge[0][1]]
                end_point = [edge[1][0], edge[1][1]]
                f["geometry"] = {"type": "LineString", "coordinates":[start_point, end_point] }
                azimuth = calculate_azimuth(start_point, end_point)
                lightness = calculate_lightness(azimuth)
                f["properties"] = {"height": height, "azimuth": azimuth, "lightness": lightness}
                mygeojson["features"].append(f)
    return mygeojson


def main():
    parser = argparse.ArgumentParser(description="My GEO1015.2024 hw01")
    parser.add_argument("inputfile", type=str, help="GeoTIFF")
    parser.add_argument(
        "thinning", type=float, help="Thinning factor (between 0 and 1)"
    )
    parser.add_argument(
        "range", type=str, help="a Python range for the contours, eg: (0, 1000, 100)"
    )

    args = parser.parse_args()

    # Validate thinning factor
    if not 0 <= args.thinning <= 1:
        parser.error("Thinning factor must be between 0 and 1")
    # Validate the range
    try:
        tmp = list(map(int, args.range.strip("() ").split(",")))
    except:
        parser.error("range invalid")
    myrange = range(tmp[0], tmp[1], tmp[2])
    print("Extracting the contours at:", list(myrange))


    try:
        # generate random points from raster
        print("Sampling points from raster to create tin...")
        sampled_points = sample_raster(args.inputfile, args.thinning)
        if not sampled_points:
            raise ValueError("Sampled points not found")
        print(f"Number of sampled points: {len(sampled_points)}")

        # create tin
        print("Creating tin...")
        dt, triangles, vertices, hull = create_tin(sampled_points)
        print("Tin created")
        print(f"Number of triangles: {len(triangles)}")

        # Generate contours
        print("\nGenerate contour lines...")
        all_intersections = get_all_intersections(dt, triangles, myrange)
        print(f"Number of intersection heights: {len(all_intersections)}")

        # organize contours
        print("\nOrganizing contour lines...")
        contours = structure_contours(all_intersections,dt)

        # create geojson
        print("Creating GeoJson...")
        geojson = create_json(contours)

        # write output file
        dt.write_ply("mydt.ply")
        print("File 'mydt.ply' created.")
        # -- write the contours to a GeoJSON file
        with open("mycontours.geojson", "w") as file:
            file.write(json.dumps(geojson, indent=2))
        print("File 'mycontour.geojson' created.")

    except Exception as e:
        print(f"Error emit : {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

