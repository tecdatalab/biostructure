import numpy as np


def get_point_face(cube, original_point, move, filter_value):
    try:
        while cube[original_point[0]][original_point[1]][original_point[2]] <= filter_value:
            original_point[0] += move[0]
            original_point[1] += move[1]
            original_point[2] += move[2]
        return original_point
    except ValueError:
        return False


def get_n_points_face(cube, points, n_points, move, filter_value, face, new_val):
    result = []
    count = 0
    actual_point = 0
    max_points = len(points)
    while count < n_points and actual_point < max_points:
        point = points[actual_point]
        if face == 0:
            point[0] = new_val
        if face == 1:
            point[1] = new_val
        if face == 2:
            point[2] = new_val
        try_point = get_point_face(cube, point, move, filter_value)
        if not isinstance(try_point, bool):
            result.append(try_point)
            count += 1
        actual_point += 1

    return result


def get_functional_points(cube, filter_value):
    pos = (np.where(cube > filter_value))

    clear_0 = np.unique(np.array([[0, pos[1][i], pos[2][i]] for i in range(len(pos[0]))]), axis=0)
    clear_1 = np.unique(np.array([[pos[0][i], 0, pos[2][i]] for i in range(len(pos[1]))]), axis=0)
    clear_2 = np.unique(np.array([[pos[0][i], pos[1][i], 0] for i in range(len(pos[2]))]), axis=0)
    return [clear_0, clear_1, clear_2]


def get_n_points_cube(cube, n_points_face, filter_value):
    values_cube = cube.shape
    points = get_functional_points(cube, filter_value)
    result = []
    np.random.shuffle(points[2])
    result += get_n_points_face(cube, points[2], n_points_face, [0, 0, 1], filter_value, 2, 0)
    np.random.shuffle(points[2])
    result += get_n_points_face(cube, points[2], n_points_face, [0, 0, -1], filter_value, 2, values_cube[2] - 1)

    np.random.shuffle(points[1])
    result += get_n_points_face(cube, points[1], n_points_face, [0, 1, 0], filter_value, 1, 0)
    np.random.shuffle(points[1])
    result += get_n_points_face(cube, points[1], n_points_face, [0, -1, 0], filter_value, 1, values_cube[1] - 1)

    np.random.shuffle(points[0])
    result += get_n_points_face(cube, points[0], n_points_face, [1, 0, 0], filter_value, 0, 0)
    np.random.shuffle(points[0])
    result += get_n_points_face(cube, points[0], n_points_face, [-1, 0, 0], filter_value, 0, values_cube[0] - 1)

    return np.unique(result, axis=0)
