# -- geo1015.2024.hw02
# -- Hugo Ledoux <h.ledoux@tudelft.nl>
# -- 2024-11-27

# ------------------------------------------------------------------------------
# DO NOT MODIFY THIS FILE!!!
# I will make my own main.py to test your `my_code_hw02.py`
# You can of course modify this file to develop the code and test things,
# but you can't write the functions taking into account a different main.
# ------------------------------------------------------------------------------

import argparse
import math
import sys

import laspy
import my_code_hw02  # -- *all* your code goes into 'my_code_hw02'
import numpy as np
import polyscope as ps
import rasterio
import startinpy


def main():
    parser = argparse.ArgumentParser(description="My GEO1015.2024 hw02")
    parser.add_argument("inputfile", type=str, help="LAZ file")
    parser.add_argument(
        "--input_thinning", type=int, default=100, help="Thinning of the input LAZ file"
    )

    parser.add_argument(
        "--resolution",
        type=int,
        default=20,
        help="GFTIN resolution for the ground grid",
    )
    parser.add_argument(
        "--max_dist", type=float, default=0.5, help="GFTIN max distance parameter"
    )
    parser.add_argument(
        "--max_angle", type=float, default=20.0, help="GFTIN max angle parameter"
    )
    args = parser.parse_args()
    # -- LAZ file exists?
    try:
        lazfile = laspy.read(args.inputfile)
    except Exception as e:
        print(e)
        sys.exit()

    # -- store non_ground points (for visualisation)
    non_ground = lazfile.classification != 2
    pts_ng = lazfile.xyz[non_ground]

    # -- bbox
    bbox = [lazfile.x.min(), lazfile.y.min(), lazfile.x.max(), lazfile.y.max()]

    print("Input LAZ file has {} points".format(lazfile.xyz.shape[0]))

    # -- DTM with classification==2
    print("=== DTM with classification ===")
    pts_ground = get_points(lazfile, classification=2, thinning=args.input_thinning)
    print("Inserting {} points".format(pts_ground.shape[0]))
    dt = startinpy.DT()
    dt.insert(pts_ground)
    grid_c_linear = rasterise("Linear", dt, bbox)
    write_grid_to_disk("dtm_c_linear.tiff", grid_c_linear, (bbox[0], bbox[3]), 1.0)
    grid_c_laplace = rasterise("Laplace", dt, bbox)
    write_grid_to_disk("dtm_c_laplace.tiff", grid_c_laplace, (bbox[0], bbox[3]), 1.0)
    # view_polyscope(dt, pts_ng)

    # -- DTM with GFTIN
    print("=== DTM with GFTIN ===")
    pts_all = get_points(lazfile, thinning=args.input_thinning)
    dt = my_code_hw02.gftin(pts_all, args.resolution, args.max_dist, args.max_angle)
    grid_gftin = rasterise("Linear", dt, bbox)
    write_grid_to_disk("dtm_gftin.tiff", grid_gftin, (bbox[0], bbox[3]), 1.0)
    # view_polyscope(dt2, pts_ng)


def get_points(lazfile, classification=None, thinning=1):
    """
    Function that reads a LAZ file and returns a numpy array Nx3 with the x, y, and z coordinates of the points.

    Inputs:
      classification:
        None: all points are used
        i:    only that class is retained
      thinning:
        1: all points are used
        10: 1/10 points are used

    Outputs:
      An Nx3 numpy array with the pointcloud (consisting of N points) that was read from the LAZ file

    """
    print("thinning factor is {}".format(thinning))
    if classification != None:
        ground = lazfile.classification == classification
        pts = lazfile.xyz[ground][::thinning]
        return pts
    else:
        pts = lazfile.xyz[::thinning]
        return pts


def rasterise(method, dt, bbox, cellsize=1.0):
    """
    Function that rasterise the DT into a (numpy) grid.

    Inputs:
      dt: the startinpy DT
      cellsize: the resolution of the grid

    Outputs:
      An NxM numpy array
    """
    deltax = math.ceil((bbox[2] - bbox[0]) / cellsize)
    deltay = math.ceil((bbox[3] - bbox[1]) / cellsize)
    centres = []
    i = 0
    for row in range((deltay - 1), -1, -1):
        j = 0
        y = bbox[1] + (row * cellsize) + (cellsize / 2)
        for col in range(deltax):
            x = bbox[0] + (col * cellsize) + (cellsize / 2)
            # -- your interpolation function is called
            if method == "Linear":
                centres.append(my_code_hw02.interpolate_linear(dt, x, y))
            elif method == "Laplace":
                centres.append(my_code_hw02.interpolate_laplace(dt, x, y))
            else:
                print("Error: interpolation method unknown")
                raise SystemExit
            j += 1
        i += 1
    centres = np.asarray(centres)
    print("Interpolated at {} locations".format(centres.shape[0]))
    return centres.reshape((deltay, deltax))


def write_grid_to_disk(output_file, a, bbox, cellsize):
    with rasterio.open(
        output_file,
        "w",
        driver="GTiff",
        height=a.shape[0],
        width=a.shape[1],
        count=1,
        dtype=np.float32,
        crs=rasterio.crs.CRS.from_string("EPSG:28992"),
        nodata=np.nan,
        transform=(cellsize, 0.0, bbox[0], 0.0, -cellsize, bbox[1]),
    ) as dst:
        dst.write(a, 1)
    print("GeoTIFF file written to '%s'" % output_file)


def view_polyscope(dt, pts_others):
    pts = dt.points
    pts[0] = pts[1]  # -- first vertex has inf and could mess things
    trs = dt.triangles
    ps.init()
    ps.set_program_name("mydt")
    ps.set_up_dir("z_up")
    ps.set_ground_plane_mode("shadow_only")
    ps.set_ground_plane_height_factor(0.01, is_relative=True)
    ps.set_autocenter_structures(True)
    ps.set_autoscale_structures(True)
    ps_mesh = ps.register_surface_mesh("mysurface", pts, trs)
    ps_mesh.reset_transform()
    pc1 = ps.register_point_cloud(
        "nonground", pts_others, radius=0.0015, point_render_mode="sphere"
    )
    pc1.reset_transform()
    ps.show()


if __name__ == "__main__":
    main()
