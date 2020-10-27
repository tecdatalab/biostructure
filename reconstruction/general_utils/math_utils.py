import math


def points_same_scale_aux(point, shape, shape_other, pos):
  if max(shape[pos], shape_other[pos]) == shape[pos]:
    return point[pos]
  else:
    return point[pos] + ((shape_other[pos] - shape[pos]) // 2)


def t_points_same_scale(point1, point2, shape1, shape2):
  new_p1 = [0, 0, 0]
  new_p2 = [0, 0, 0]

  new_p1[0] = points_same_scale_aux(point1, shape1, shape2, 0)
  new_p1[1] = points_same_scale_aux(point1, shape1, shape2, 1)
  new_p1[2] = points_same_scale_aux(point1, shape1, shape2, 2)

  new_p2[0] = points_same_scale_aux(point2, shape2, shape1, 0)
  new_p2[1] = points_same_scale_aux(point2, shape2, shape1, 1)
  new_p2[2] = points_same_scale_aux(point2, shape2, shape1, 2)

  return new_p1, new_p2


def get_vector_move_1to2(point1, point2):
  result = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
  return result


def chance_base_point(point, actual_shape, new_shape):
  result = [0, 0, 0]

  result[0] = new_shape[0] * (point[0] / actual_shape[0])
  result[1] = new_shape[1] * (point[1] / actual_shape[1])
  result[2] = new_shape[2] * (point[2] / actual_shape[2])

  return result


def distance_3d_points(point1, point2):
  dist = math.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2 + (point1[2] - point2[2]) ** 2)
  return dist
