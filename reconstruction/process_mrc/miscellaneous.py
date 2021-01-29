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


