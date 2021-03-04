import numpy as np


def get_center_point(index_segments, segments, filter_val):
  points_pos = None

  for i in index_segments:
    for j in segments:
      if j.id_segment == i:
        actual_segment = j
        temp_points = np.where(actual_segment.mask > filter_val)
        if points_pos is None:
          points_pos = np.column_stack(temp_points)
        else:
          points_pos = np.concatenate((points_pos, np.column_stack(temp_points)), axis=0)

  center_point = np.mean(points_pos, axis=0)
  return [int(center_point[0]), int(center_point[1]), int(center_point[2])]


def get_center_point_by_graph(index_segments, graph):
  x = 0
  y = 0
  z = 0
  total_div_points = 0

  for i in index_segments:
    x += graph.nodes[i]["cube_xyz_can"][0][0]
    y += graph.nodes[i]["cube_xyz_can"][0][1]
    z += graph.nodes[i]["cube_xyz_can"][0][2]
    total_div_points += graph.nodes[i]["cube_xyz_can"][1]

  return [int(x / total_div_points), int(y / total_div_points), int(z / total_div_points)]


def get_sum_xyz_can(segment):
  points = np.where(segment.mask > 0)
  points_pos = np.column_stack(points)
  total_sum = np.sum(points_pos, axis=0)

  return [total_sum.tolist(), points_pos.size]
