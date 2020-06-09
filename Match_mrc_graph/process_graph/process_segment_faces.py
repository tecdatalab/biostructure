from random import seed
from random import randint
import numpy as np

def get_point_face(cube, original_point, move, filter_value):
    try:
        while cube[original_point[0]][original_point[1]][original_point[2]]>filter_value:
            original_point[0]+= move[0]
            original_point[1]+= move[1]
            original_point[2]+= move[2]
        return original_point
    except:
        return False

def get_n_points_face(cube, points, n_points, move, filter_value):
    result = []
    count = 0
    actual_point = 0
    max_points =len(points)
    while count < n_points and actual_point<max_points:
        point=points[actual_point]
        try_point = get_point_face(cube, point, move, filter_value)
        if not isinstance(try_point, bool):
            result.append(try_point)
            count+=1
        actual_point+=1
            
    return result

def get_funtional_points(pos, cube, filter_value, new_value):
    pos = (np.where(cube > filter_value))
    clear = None
    if pos==0:
        clear = np.unique(np.array([[new_value,pos[1][i],pos[2][i]] for i in range(len(pos[0]))]), axis=0)
    elif pos==1:
        clear = np.unique(np.array([[pos[0][i],new_value,pos[2][i]] for i in range(len(pos[0]))]), axis=0)
    else:
        clear = np.unique(np.array([[pos[0][i],pos[1][i],new_value] for i in range(len(pos[0]))]), axis=0)
    np.random.shuffle(clear)
    return clear

def get_n_points_cube(cube, n_points_face, filter_value):
    values_cube = cube.shape
    result = []
    
    points = get_funtional_points(2, cube, filter_value, 0)
    result+= get_n_points_face(cube,points, n_points_face, [0,0,1],filter_value)
    points = get_funtional_points(2, cube, filter_value, values_cube[2]-1)
    result+= get_n_points_face(cube,points, n_points_face, [0,0,-1],filter_value)
    
    points = get_funtional_points(1, cube, filter_value, 0)
    result+= get_n_points_face(cube,points, n_points_face, [0,1,0], filter_value)
    points = get_funtional_points(1, cube, filter_value, values_cube[1]-1)
    result+= get_n_points_face(cube,points, n_points_face, [0,-1,0], filter_value)
    
    points = get_funtional_points(0, cube, filter_value, 0)
    result+= get_n_points_face(cube,points, n_points_face, [1,0,0], filter_value)
    points = get_funtional_points(0, cube, filter_value, values_cube[0]-1)
    result+= get_n_points_face(cube,points, n_points_face, [-1,0,0], filter_value)
    
    return result
       
