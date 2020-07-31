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
    return (int(np.median(x_pos)), int(np.median(y_pos)), int(np.median(z_pos)))
