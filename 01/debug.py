import numpy as np
import rasterio
import startinpy
import json
from geo1015_hw01 import sample_raster,create_tin,get_intersection,get_all_intersections,create_json,structure_contours
# 引用你的函数代码
# 确保你将函数定义写在一个文件中并导入它们，例如 `from my_functions import sample_raster, create_tin, get_all_intersections`

# 测试代码
def main():
    # 输入参数
    dataset = "dem_01.tif"  # 替换为实际的 GeoTIFF 文件路径
    thinning = 0.01  # 采样比例
    heights = range(100,1000,100)  # 等高线高度值列表

    print("Step 1: Sampling raster data...")
    sampled_points = sample_raster(dataset, thinning)
    print(f"Sampled Points: {len(sampled_points)} points")

    print("Step 2: Creating TIN...")
    dt, triangles, vertices, hull = create_tin(sampled_points)
    print(f"TIN created with {len(triangles)} triangles and {len(vertices)} vertices")

    print("Step 3: Calculating intersections...")
    intersections = get_all_intersections(dt, triangles, heights)
    # with open('intersections.txt', "w") as f:
    #     f.write(str(intersections))
    print("Step 4: Organizing intersections...")
    organized_contours = structure_contours(intersections,dt)
    result = create_json(organized_contours)
    with open("test.geojson", "w") as file:
        file.write(json.dumps(result, indent=2))
    print("File 'test.geojson' created.")
    dt.write_ply("testdt.ply")
    print("File 'testdt.ply' created.")



if __name__ == "__main__":
    main()

