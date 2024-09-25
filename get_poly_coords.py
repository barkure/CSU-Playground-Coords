import math
import numpy as np
import matplotlib.pyplot as plt
from pygcj.pygcj import GCJProj


def calc_coordinates_by_deg(lon, lat, deg, distance):
    """
    Calculate new coordinates based on the given distance and degree
    :param lon: longitude
    :param lat: latitude
    :param deg: degree
    :param distance: distance
    :return: new coordinates
    """
    earth_radius = 6371.393 * 1000 # 地球半径
    meters_per_degree = earth_radius * math.pi / 180 # 每度对应的米数（赤道上）
    rad = math.radians(deg)
    new_lon = lon + (distance * math.sin(rad)) / (meters_per_degree * math.cos(math.radians(lat)))
    new_lat = lat + (distance * math.cos(rad)) / meters_per_degree
    return new_lon, new_lat


def calc_coordinates_by_xy(lon, lat, x, y):
    """
    Calculate new coordinates based on the given x and y
    :param lon: longitude
    :param lat: latitude
    :param x: x
    :param y: y
    :return: new coordinates
    """
    earth_radius = 6371.393 * 1000 # 地球半径
    meters_per_degree = earth_radius * math.pi / 180 # 每度对应的米数（赤道上）
    new_lon = lon + x / (meters_per_degree * math.cos(math.radians(lat)))
    new_lat = lat + y / meters_per_degree
    return new_lon, new_lat


def rotate_point(point, angle, center):
    """
    Rotate a point around a center
    :param point: point to rotate
    :param angle: angle to rotate in radians
    :param center: center of rotation
    :return: rotated point
    """
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    rotation_matrix = np.array([
        [cos_theta, -sin_theta],
        [sin_theta, cos_theta]
    ])
    translated_point = point - center
    rotated_point = np.dot(rotation_matrix, translated_point)
    final_point = rotated_point + center
    return final_point


# 创建操场多边形（中南大学新校区体育场副田径场）
"""
以上侧（北方）的弧（此处当作半圆）的圆心为基准点，该点坐标为 (112.927251362, 28.160720451)
基准点南偏东 22°（方位角为 158°），距离 79 米的点为下侧半圆的圆心,半圆半径均为 40 米

然后凭感觉用直线连接两个半圆的端点，形成一个操场多边形

绘制倾斜的图形对猪脑要求太大了，我们先绘制图形，然后离散化点，最后将点旋转到指定方位就好了
"""

# 定义半圆的圆心和半径
center1 = (0, 0)
center2 = (0, -79)
radius = 40

# 计算半圆的弧长和总的周长
arc_length = math.pi * radius
total_length = arc_length * 2 + 79 * 2

# 计算每个点的间距
num_points = 100
spacing = total_length / num_points

# 生成等距离的点
points = []
current_distance = 0

# 生成上半圆上的点
for i in range(num_points):
    if current_distance < arc_length:
        angle = current_distance / radius
        x = center1[0] + radius * math.cos(angle)
        y = center1[1] + radius * math.sin(angle)
        points.append([x, y])
    else:
        break
    current_distance += spacing

# 生成左侧直线上的点
line_start = (center1[0] - radius, center1[1])
line_end = (center2[0] - radius, center2[1])
for i in range(num_points):
    if current_distance < arc_length + 79:
        t = (current_distance - arc_length) / 79
        x = line_start[0] + t * (line_end[0] - line_start[0])
        y = line_start[1] + t * (line_end[1] - line_start[1])
        points.append([x, y])
    else:
        break
    current_distance += spacing

# 生成下半圆上的点
for i in range(num_points):
    if current_distance < arc_length * 2 + 79:
        angle = (current_distance - arc_length - 79) / radius
        x = center2[0] + radius * math.cos(angle + math.pi)
        y = center2[1] + radius * math.sin(angle + math.pi)
        points.append([x, y])
    else:
        break
    current_distance += spacing

# 生成右侧直线上的点
line_start = (center2[0] + radius, center2[1])
line_end = (center1[0] + radius, center1[1])
for i in range(num_points):
    if current_distance < total_length:
        t = (current_distance - arc_length * 2 - 79) / 79
        x = line_start[0] + t * (line_end[0] - line_start[0])
        y = line_start[1] + t * (line_end[1] - line_start[1])
        points.append([x, y])
    else:
        break
    current_distance += spacing

points = np.array(points)

# 以上侧半圆圆心为基准点，旋转 22°
angle = math.radians(22)
rotate_points = np.array([rotate_point(point, angle, center1) for point in points])

# 绘制闭合图形和离散点
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(points[:, 0], points[:, 1], 'b-', linewidth=1, label='Discrete Points')
plt.scatter(points[:, 0], points[:, 1], color='red', s=5, label='Discrete Points')
plt.axis('equal')
plt.grid(True)
plt.xlabel('X (m)', fontsize=14, fontweight='bold')
plt.ylabel('Y (m)', fontsize=14, fontweight='bold')

# 绘制旋转后的图形
plt.subplot(1, 2, 2)
plt.plot(rotate_points[:, 0], rotate_points[:, 1], 'b-', linewidth=1, label='Discrete Points')
plt.scatter(rotate_points[:, 0], rotate_points[:, 1], color='red', s=5, label='Discrete Points')
plt.axis('equal')
plt.grid(True)
plt.xlabel('X (m)', fontsize=14, fontweight='bold')
plt.ylabel('Y (m)', fontsize=14, fontweight='bold')

plt.savefig('points.png')


# 基准点经纬度坐标
coord1 = (112.927251362, 28.160720451)
# 将点转换为经纬度坐标
coords = np.array([calc_coordinates_by_xy(coord1[0], coord1[1], point[0], point[1]) for point in rotate_points])
# 导出为 CSV 文件
np.savetxt('coords.csv', coords, delimiter=',', fmt='%.10f')


# 将坐标转换为 GCJ-02 坐标
trans = GCJProj()
gcj_coords = np.array([trans.wgs_to_gcj(coord[1], coord[0]) for coord in coords])
# 导出为 CSV 文件（交换两列）
np.savetxt('coords_gcj.csv', gcj_coords[:, [1, 0]], delimiter=',', fmt='%.10f')
