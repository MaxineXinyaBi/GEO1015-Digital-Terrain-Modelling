def get_all_intersections(dt,triangles,heights,tolerance=1e-10):
    intersections_by_height = {height: {"closed": [], "open": []} for height in heights}
    points = dt.points[1:]
    for triangle in triangles:
        v1 = points[triangle[0] - 1]
        v2 = points[triangle[1] - 1]
        v3 = points[triangle[2] - 1]

        z_min = min(v1[2], v2[2], v3[2]) - tolerance
        z_max = max(v1[2], v2[2], v3[2]) + tolerance

        for height in heights:
            if height < z_min or height > z_max:
                continue

            # 计算三角形与等高线的交点
            intersections = get_intersection(dt, triangle, height, tolerance=tolerance)

            # 如果存在交点，将其添加到对应高程的列表中
            if intersections:
                # 判断是否为闭合线段
                if intersections[0] == intersections[-1]:  # 闭合线段
                    intersections_by_height[height]["closed"].append(intersections)
                else:  # 开放线段
                    intersections_by_height[height]["open"].append(intersections)
    # 对交点进行方向统一处理
    for height in intersections_by_height:
        contours = intersections_by_height[height]
        # 处理闭合等高线方向
        contours["closed"] = [ensure_counter_clockwise(c) for c in contours["closed"]]
        # 处理开放等高线方向
        contours["open"] = [ensure_open_direction(c) for c in contours["open"]]

    return intersections_by_height


def is_counter_clockwise(points):
    """
    判断闭合线段是否是逆时针方向。

    Args:
        points: 点列表 [(x, y), ...]，闭合线段第一个点和最后一个点相同。

    Returns:
        bool: True 如果是逆时针方向，False 如果是顺时针方向。
    """
    n = len(points)
    area = 0
    for i in range(n - 1):  # 不包括最后一个点，最后一个点与第一个点计算
        x1, y1 = points[i]
        x2, y2 = points[i + 1]
        area += (x2 - x1) * (y2 + y1)
    return area > 0  # 如果面积为正，则是逆时针方向


def ensure_counter_clockwise(points):
    """
    确保闭合线段是逆时针方向。

    Args:
        points: 点列表 [(x, y), ...]，闭合线段第一个点和最后一个点相同。

    Returns:
        list: 如果是顺时针方向，则反转点顺序以调整为逆时针。
    """
    if not is_counter_clockwise(points):
        points.reverse()  # 如果是顺时针方向，反转点顺序
    return points


def ensure_open_direction(points):
    """
    确保开放线段的方向一致。

    Args:
        points: 点列表 [(x, y), ...]，开放线段。

    Returns:
        list: 调整方向后的一致点列表。
    """
    start, end = points[0], points[-1]
    if start > end:  # 假设起点需要小于终点（按照 x 或 y 排序）
        points.reverse()  # 如果方向不一致，则反转点顺序
    return points