import numpy as np
import startinpy
import my_code_hw02
import laspy
import random
import time


def load_point_cloud(file_path, classification=None, thinning=100):
    """
    加载点云文件并返回选定的点

    Parameters:
    file_path: 点云文件路径
    classification: 可选，指定点的分类（如地面点classification=2）
    thinning: 抽稀因子，每隔多少个点取一个
    """
    try:
        las = laspy.read(file_path)
        if classification is not None:
            # 只选择特定分类的点
            mask = las.classification == classification
            points = las.xyz[mask][::thinning]
        else:
            # 选择所有点
            points = las.xyz[::thinning]

        return points
    except Exception as e:
        print(f"Error loading point cloud file: {e}")
        return None


def create_test_points(bbox, num_points=10):
    """
    在边界框内随机生成测试点
    """
    x_min, y_min = bbox[0], bbox[1]
    x_max, y_max = bbox[2], bbox[3]

    test_points = []
    for _ in range(num_points):
        x = random.uniform(x_min, x_max)
        y = random.uniform(y_min, y_max)
        test_points.append([x, y])

    return test_points


def compare_interpolations(dt, test_points):
    """
    比较不同插值方法的结果
    """
    results = []

    for i, (x, y) in enumerate(test_points):
        print(f"\nTest point {i + 1}: ({x:.2f}, {y:.2f})")

        # 我们的实现
        start_time = time.time()
        our_linear = my_code_hw02.interpolate_linear(dt, x, y)
        linear_time = time.time() - start_time

        start_time = time.time()
        our_laplace = my_code_hw02.interpolate_laplace(dt, x, y)
        laplace_time = time.time() - start_time

        # startinpy的实现（仅用于比较）
        start_time = time.time()
        startinpy_linear = dt.interpolate({"method": "TIN"}, [[x, y]])[0]
        startinpy_linear_time = time.time() - start_time

        start_time = time.time()
        startinpy_laplace = dt.interpolate({"method": "Laplace"}, [[x, y]])[0]
        startinpy_laplace_time = time.time() - start_time

        results.append({
            'point': (x, y),
            'our_linear': our_linear,
            'our_laplace': our_laplace,
            'startinpy_linear': startinpy_linear,
            'startinpy_laplace': startinpy_laplace,
            'linear_diff': abs(our_linear - startinpy_linear),
            'laplace_diff': abs(our_laplace - startinpy_laplace),
            'our_linear_time': linear_time,
            'our_laplace_time': laplace_time,
            'startinpy_linear_time': startinpy_linear_time,
            'startinpy_laplace_time': startinpy_laplace_time
        })

        print(f"Linear interpolation:")
        print(f"  Our implementation: {our_linear:.3f} (Time: {linear_time:.4f}s)")
        print(f"  Startinpy: {startinpy_linear:.3f} (Time: {startinpy_linear_time:.4f}s)")
        print(f"  Difference: {abs(our_linear - startinpy_linear):.6f}")

        print(f"Laplace interpolation:")
        print(f"  Our implementation: {our_laplace:.3f} (Time: {laplace_time:.4f}s)")
        print(f"  Startinpy: {startinpy_laplace:.3f} (Time: {startinpy_laplace_time:.4f}s)")
        print(f"  Difference: {abs(our_laplace - startinpy_laplace):.6f}")

    return results


def analyze_results(results):
    """
    分析测试结果
    """
    linear_diffs = [r['linear_diff'] for r in results]
    laplace_diffs = [r['laplace_diff'] for r in results]

    print("\nAnalysis Summary:")
    print(f"Linear interpolation differences:")
    print(f"  Max: {max(linear_diffs):.6f}")
    print(f"  Mean: {np.mean(linear_diffs):.6f}")
    print(f"  Std: {np.std(linear_diffs):.6f}")

    print(f"Laplace interpolation differences:")
    print(f"  Max: {max(laplace_diffs):.6f}")
    print(f"  Mean: {np.mean(laplace_diffs):.6f}")
    print(f"  Std: {np.std(laplace_diffs):.6f}")


def main():
    # 设置随机种子以保证可重复性
    random.seed(42)

    # 加载点云文件
    point_cloud_file = "crop.laz"  # 替换为你的点云文件路径
    points = load_point_cloud(point_cloud_file, classification=2, thinning=100)

    if points is None:
        return

    print(f"Loaded {len(points)} points from point cloud")

    # 创建Delaunay三角网
    dt = startinpy.DT()
    dt.insert(points)

    # 获取边界框
    bbox = [
        points[:, 0].min(), points[:, 1].min(),  # x_min, y_min
        points[:, 0].max(), points[:, 1].max()  # x_max, y_max
    ]

    # 生成测试点
    test_points = create_test_points(bbox, num_points=10)

    # 运行测试
    print("\nRunning interpolation tests...")
    results = compare_interpolations(dt, test_points)

    # 分析结果
    analyze_results(results)


if __name__ == "__main__":
    main()