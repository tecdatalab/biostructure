import numpy as np


def get_center_point(index_segments, segments, filter):
    x_pos = np.array([])
    y_pos = np.array([])
    z_pos = np.array([])

    con = 0
    for i in index_segments:
        for j in segments:
            if j.id_segment == i:
                actual_segment = j
                temp_points = np.where(actual_segment.mask > filter)
                x_pos = np.concatenate((x_pos, temp_points[0]), axis=0)
                y_pos = np.concatenate((y_pos, temp_points[1]), axis=0)
                z_pos = np.concatenate((z_pos, temp_points[2]), axis=0)
                break
    return [int(np.median(x_pos)), int(np.median(y_pos)), int(np.median(z_pos))]


def transform_points_sscale_aux(point, shape, shape_other, pos):
    if max(shape[pos], shape_other[pos]) == shape[pos]:
        return point[pos]
    else:
        return point[pos] + ((shape_other[pos] - shape[pos]) // 2)


def transform_points_sscale(point1, point2, shape1, shape2):
    newp1 = [0, 0, 0]
    newp2 = [0, 0, 0]

    newp1[0] = transform_points_sscale_aux(point1, shape1, shape2, 0)
    newp1[1] = transform_points_sscale_aux(point1, shape1, shape2, 1)
    newp1[2] = transform_points_sscale_aux(point1, shape1, shape2, 2)

    newp2[0] = transform_points_sscale_aux(point2, shape2, shape1, 0)
    newp2[1] = transform_points_sscale_aux(point2, shape2, shape1, 1)
    newp2[2] = transform_points_sscale_aux(point2, shape2, shape1, 2)

    return newp1, newp2


def get_vector_move_1to2(point1, point2):
    result = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
    return result
