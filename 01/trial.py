import argparse
import json
import sys
import rasterio
import startinpy
import numpy as np
import math

def main():
    parser = argparse.ArgumentParser(description="My GEO1015.2024 hw01")
    parser.add_argument("inputfile", type=str, help="GeoTIFF")
    parser.add_argument(
        "thinning", type=float, help="Thinning factor (between 0 and 1)"
    )
    parser.add_argument(
        "range", type=str, help="a Python range for the contours, eg: (100, 500, 100)"
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

    # -- load in memory the input GeoTIFF
    try:
        d = rasterio.open(args.inputfile)
    except Exception as e:
        print(e)
        sys.exit()

    print("name:", d.name)
    print("crs:", d.crs)
    print("size:", d.shape)
    print("no_data:", d.nodata)

    dtm = d.read(1)  # Read elevation data
    sampled_points = sample_dtm(dtm, args.thinning, d.transform)

    # Create TIN
    dt = startinpy.DT()
    for point in sampled_points:
        dt.insert_one_pt(point)
    print(f"Inserted {len(sampled_points)} points into the TIN.")

    # Visit triangles
    visit_triangles(dt)

    # Extract contours
    contours = extract_contours(dt, list(myrange))
    contours = add_tanaka_properties(contours)
    contours = enforce_counter_clockwise_orientation(contours)
    print(f"Extracted {len(contours)} contours.")

    # -- write the triangulation to a PLY file
    dt.write_ply("mydt.ply")
    print("File 'mydt.ply' created.")

    # -- write the contours to a GeoJSON file
    write_geojson(contours, "mycontours.geojson")
    print("File 'mycontours.geojson' created.")


def sample_dtm(dtm, thinning, transform):
    rows, cols = dtm.shape
    sampled_points = []
    for r in range(rows):
        for c in range(cols):
            if np.random.random() < thinning:  # Probabilistically sample points
                z = dtm[r, c]
                if z != -9999:  # Avoid no-data values
                    x, y = rasterio.transform.xy(transform, r, c)
                    sampled_points.append([x, y, z])
    return sampled_points


def extract_contours(dt, contour_range):
    contours = []
    trs = dt.triangles
    for z in contour_range:
        lines = dt.extract_contours(z)  # Extract isolines at elevation z
        for line in lines:
            contours.append({"elevation": z, "coordinates": line})

    return contours

def visit_triangles(dt):
    trs = dt.triangles
    triangles_to_remove = []
    for i, triangle in enumerate(trs):
        # Get the coordinates of the triangle vertices
        v1, v2, v3 = triangle
        p1 = dt.points[v1]
        p2 = dt.points[v2]
        p3 = dt.points[v3]

        # Check if the triangle is horizontal
        if p1[2] == p2[2] == p3[2]:
            # If the triangle is horizontal, mark it for removal
            triangles_to_remove.append(i)
        else:
            pass

    # Remove horizontal triangles
    for index in sorted(triangles_to_remove, reverse=True):
        del dt.triangles[index]

def calculate_azimuth(coordinates):
    x1, y1 = coordinates[0]
    x2, y2 = coordinates[-1]
    angle = math.atan2(y2 - y1, x2 - x1) * 180 / math.pi
    return (angle + 360) % 360  # Normalize to [0, 360]

def calculate_lightness(azimuth, light_source=(315, 45)):  # NW by default
    l = abs(((azimuth - 225) % 360) - 180)
    lightness = (l / 180) * 100
    return round(lightness, 2)

def add_tanaka_properties(contours, light_source=(315, 45)):
    for contour in contours:
        azimuth = calculate_azimuth(contour["coordinates"])
        lightness = calculate_lightness(azimuth, light_source)
        contour.update({"azimuth": azimuth, "lightness": lightness})
    return contours


def enforce_counter_clockwise_orientation(contours):
    for contour in contours:
        if is_clockwise(contour["coordinates"]):
            contour["coordinates"].reverse()  # Reverse the order of coordinates
    return contours


def is_clockwise(coords):
    total = 0
    for i in range(len(coords)):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % len(coords)]
        total += (x2 - x1) * (y2 + y1)
    return total > 0  # Clockwise if positive


def write_geojson(contours, filename):
    geojson = {"type": "FeatureCollection", "features": []}
    for contour in contours:
        feature = {
            "type": "Feature",
            "geometry": {"type": "LineString", "coordinates": contour["coordinates"]},
            "properties": {
                "elevation": contour["elevation"],
                "azimuth": contour.get("azimuth"),
                "lightness": contour.get("lightness"),
            },
        }
        geojson["features"].append(feature)

    with open(filename, "w") as file:
        json.dump(geojson, file, indent=2)
        print(f"GeoJSON saved to '{filename}'.")

    return geojson

if __name__ == "__main__":
    main()