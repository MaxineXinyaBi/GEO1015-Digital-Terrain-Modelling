def create_json(contours):
    geojson = {
        "type": "FeatureCollection",
        "features": []
    }

    for height, contour_groups in contours.items():
        for type_key in ["closed", "open"]:
            for contour in contour_groups[type_key]:
                for segment in contour:
                    azimuth = calculate_azimuth(segment[0], segment[1])
                    lightness = calculate_lightness(azimuth)

                    feature = {
                        "type": "Feature",
                        "geometry": {
                            "type": "LineString",
                            "coordinates": [
                                [segment[0][0], segment[0][1]],
                                [segment[1][0], segment[1][1]]
                            ]
                        },
                        "properties": {
                            "height": height,
                            "azimuth": azimuth,
                            "lightness": lightness
                        }
                    }
                    geojson["features"].append(feature)

    return geojson


# modified geojson
# def create_json(all_intersections):
#     geojson = {
#         "type": "FeatureCollection",
#         "features": []
#     }
#
#     total_segments = sum(len(segments) for segments in all_intersections.values())
#     print(f"Processing {total_segments} segments across {len(all_intersections)} height levels")
#
#     for height, segments in all_intersections.items():
#         for segment in segments:
#             azimuth = calculate_azimuth(segment[0], segment[1])
#             lightness = calculate_lightness(azimuth)
#
#             feature = {
#                 "type": "Feature",
#                 "geometry": {
#                     "type": "LineString",
#                     "coordinates": [
#                         [segment[0][0], segment[0][1]],
#                         [segment[1][0], segment[1][1]]
#                     ]
#                 },
#                 "properties": {
#                     "height": float(height),
#                     "azimuth": azimuth,
#                     "lightness": lightness
#                 }
#             }
#             geojson["features"].append(feature)
#
#     print(f"Created GeoJSON with {len(geojson['features'])} features")
#     return geojson


# def organize_contours(all_intersections, tolerance=1e-10):
#     """organize contours line segment on different heights based on ccw order"""
#     organized_contours = {}
#
#     # 测试用
#     total_heights = len(all_intersections)
#     print(f"\nOrganizing contours for {total_heights} height levels...")
#     for i, (height, segments) in enumerate(all_intersections.items()):
#         print(f"Processing height {height} ({i + 1}/{total_heights})")
#         print(f"Number of segments at this height: {len(segments)}")
#
#
#     for height, segments in all_intersections.items():
#         remaining_segments = segments[:]
#         organized_contours[height] = {"closed":[],"open":[]}
#
#         while remaining_segments:
#             # 测试用
#             if len(remaining_segments) % 100 == 0:  # 每处理100个段打印一次进度
#                 print(f"  Remaining segments: {len(remaining_segments)}")
#
#
#             current_segment = remaining_segments.pop(0)
#             # current_contour = [[current_segment[0], current_segment[1]], ]
#             current_segment = [(float(p[0]), float(p[1])) for p in current_segment]
#             current_contour = [current_segment]
#
#             while True:
#                 connected = False
#                 for i, segment in enumerate(remaining_segments):
#                     start, end = segment
#                     if check_connected_points(current_contour[-1][-1], start, tolerance):
#                         current_contour.append([start, end])
#                         remaining_segments.pop(i)
#                         connected = True
#                         break
#                 if not connected:
#                     break
#
#             if check_connected_points(current_contour[0][0], current_contour[-1][1], tolerance):
#                 organized_contours[height]["closed"].append(current_contour)
#             else:
#                 organized_contours[height]["open"].append(current_contour)
#
#         # 测试用
#         print(f"Finished height {height}: {len(organized_contours[height]['closed'])} closed contours, "
#               f"{len(organized_contours[height]['open'])} open contours")
#
#     return organized_contours



def get_all_intersections(dt, triangles, heights, tolerance=1e-10):
    """organize all the intersections into a dict"""
    all_intersections = {height:[] for height in heights}

    # 测试用
    total_triangles = len(triangles)
    print(f"Processing {total_triangles} triangles for {len(heights)} height levels...")
    points = dt.points[1:]

    for i, triangle in enumerate(triangles):
        if i % 100 == 0:  # 每处理1000个三角形打印一次进度
            print(f"Progress: {i}/{total_triangles} triangles processed ({(i / total_triangles * 100):.1f}%)")

        triangle_points = [points[i-1] for i in triangle]
        min_z = min(p[2] for p in triangle_points)
        max_z = max(p[2] for p in triangle_points)
        possible_heights = [h for h in heights if min_z-tolerance <= h <= max_z+tolerance]

        for height in possible_heights:
            intersection_points = get_intersection(dt, triangle, height, tolerance)
            if intersection_points:
                if i % 1000 == 0:
                    print(f"Found intersection at height {height}: {points}")
                for i in range(0, len(intersection_points), 2):
                    if i + 1 < len(intersection_points):
                        segment = [
                            (intersection_points[i][0], intersection_points[i][1]),
                            (intersection_points[i + 1][0], intersection_points[i + 1][1])
                        ]
                    all_intersections[height].append(segment)

    # 打印结果统计
    print("\nIntersection statistics:")
    for height in heights:
        print(f"Height {height}: {len(all_intersections[height])} intersections")

    return all_intersections

def check_connected_points(p1, p2, tolerance):
    """ check if two points are connected to each other within a tolerance"""
    if abs(p1[0] - p2[0]) <= tolerance and abs(p1[1] - p2[1]) <= tolerance:
        return True
    else:
        return False

def organize_contours(all_intersections, tolerance=1e-10):
    """organize contours line segment on different heights based on ccw order"""
    organized_contours = {}

    for height, segments in all_intersections.items():
        remaining_segments = segments[:]  # 复制一份避免修改原始数据
        organized_contours[height] = {"closed": [], "open": []}

        while remaining_segments:
            current_segment = remaining_segments.pop(0)
            current_contour = [current_segment]

            searching = True
            while searching:
                searching = False
                i = 0
                while i < len(remaining_segments):
                    segment = remaining_segments[i]
                    next_points = [(float(p[0]), float(p[1])) for p in segment]

                    # 检查当前轮廓末端是否能与下一段相连
                    if check_connected_points(current_contour[-1][-1], next_points[0], tolerance):
                        current_contour.append(next_points)
                        remaining_segments.pop(i)
                        searching = True
                        break
                    # 检查是否需要反转下一段
                    elif check_connected_points(current_contour[-1][-1], next_points[-1], tolerance):
                        next_points.reverse()
                        current_contour.append(next_points)
                        remaining_segments.pop(i)
                        searching = True
                        break
                    i += 1

            # 检查是否形成闭合轮廓
            if check_connected_points(current_contour[0][0], current_contour[-1][-1], tolerance):
                organized_contours[height]["closed"].append(current_contour)
            else:
                organized_contours[height]["open"].append(current_contour)

    return organized_contours
