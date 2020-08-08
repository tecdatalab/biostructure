'''
Created on 26 jul. 2020

@author: dnnxl
'''
from psycopg2 import sql

class Segment_entry(object):
    '''
    The Segment_entry object to store the segment information.

    Attributes:
        __segment_id : int
            Store the segment id.
        __segment_map_id : int
            Store the map id related to the segment.
        __segment_id_algorithm : int
            Store the id of the algorithm.
        __segment_volume : float
            Store the volume of the segment.
        __segment_contour_level : int
            Store the contour level of the segment.
        __segment_path : string
            Store the path of the segment map generate.
    '''
    __segment_id = None
    __segment_map_id = None
    __segment_id_algorithm = None
    __segment_volume = None
    __segment_contour_level = None
    __segment_path = None

    def __init__(
            self, 
            segment_id, 
            segment_map_id, 
            segment_id_algorithm,
            segment_volume,
            segment_contour_level,
            segment_path):
        self.__segment_id = segment_id
        self.__segment_map_id = segment_map_id
        self.__segment_id_algorithm = segment_id_algorithm
        self.__segment_volume = segment_volume
        self.__segment_contour_level = segment_contour_level
        self.__segment_path = segment_path
    
    def get_segment_id(self):
        return self.__segment_id
    
    def get_segment_map_id(self):
        return self.__segment_map_id
    
    def get_segment_id_algorithm(self):
        return self.__segment_id_algorithm
    
    def get_segment_volume(self):
        return self.__segment_volume

    def get_segment_contour_level(self):
        return self.__segment_contour_level
    
    def get_segment_path(self):
        return self.__segment_path

    def set_segment_id(self, value):
        self.__segment_id = value
    
    def set_segment_map_id(self, value):
        self.__segment_map_id = value
    
    def set_segment_id_algorithm(self, value):
        self.__segment_id_algorithm = value
    
    def set_segment_volume(self, value):
        self.__segment_volume = value

    def set_segment_contour_level(self, value):
        self.__segment_contour_level = value
    
    def set_segment_path(self, value):
        self.__segment_path = value

    def del_segment_id(self):
        del self.__segment_id 
    
    def del_segment_map_id(self):
        del self.__segment_map_id 
    
    def del_segment_id_algorithm(self):
        del self.__segment_id_algorithm
    
    def del_segment_volume(self):
        del self.__segment_volume 

    def del_segment_contour_level(self):
        del self.__segment_contour_level 

    def del_segment_path(self):
        del self.__segment_path 

    segment_id = property(
        get_segment_id, 
        set_segment_id, 
        del_segment_id, 
        "segment_id's docstring")
    segment_map_id = property(
        get_segment_map_id,
        set_segment_map_id,
        del_segment_map_id,
        "segment_map_id's docstring")
    segment_id_algorithm = property(
        get_segment_id_algorithm,
        set_segment_id_algorithm,
        del_segment_id_algorithm,
        "segment_id_algorithm's docstring")
    segment_volume = property(
        get_segment_volume,
        set_segment_volume,
        del_segment_volume,
        "segment_volume's docstring")
    segment_contour_level = property(
        get_segment_contour_level,
        set_segment_contour_level,
        del_segment_contour_level,
        "segment_contour_level's docstring")
    segment_path = property(
        get_segment_path,
        set_segment_path,
        del_segment_path,
        "segment_path's docstring")

    def __insert_segment_entry(self, cur):
        cur.execute(
            sql.SQL("INSERT INTO segment_entry(map_id, algorithm_id, countour_level, volume, path_map) VALUES(%s, %s, %s, %s, %s) RETURNING id_segment;"), 
                [self.__segment_map_id, self.__segment_id_algorithm, self.__segment_contour_level, self.__segment_volume, self.__segment_path])
        self.__segment_id = [record for record in cur][0]
    
    def __update_segment_entry(self, cur):
        cur.execute(sql.SQL("UPDATE segment_entry SET map_id = %s, algorithm_id = %s, countour_level = %s, volume = %s, path_map = %s WHERE id_segment = %s;"), [
            self.__segment_map_id, self.__segment_id_algorithm, self.__segment_contour_level, self.__segment_volume, self.__segment_path, self.__segment_id_algorithm])

    def insert_db(self, cur):
        self.__insert_segment_entry(cur)

    def update_db(self, cur):
        self.__update_segment_entry(cur)

    def insert_update_db(self, cur):
        cur.execute(
            sql.SQL("SELECT id_segment FROM segment_entry WHERE map_id = %s AND algorithm_id = %s AND countour_level = %s AND volume = %s"), [
                self.__segment_map_id, self.__segment_id_algorithm, self.__segment_contour_level, self.__segment_volume])
        result = [record[0] for record in cur]
        if len(result)>0:
            self.set_segment_id(result[0])
            self.update_db(cur)
        else:
            self.insert_db(cur)

class List_segments(object):
    '''
    The List_segments object to store the list of segment information.

    Attributes:
        __list_segments : list of Segment_entry object.
            Store the list of segments object.
    '''
    __list_segments = None

    def __init__(
            self, list_segments):
        self.__list_segments = list_segments

    def get_list_segments(self):
        return self.__list_segments
    
    def set_list_segments(self, value):
        self.__list_segments = value
    
    def del_list_segments(self):
        del self.__list_segments
    
    list_segments = property(
        get_list_segments, 
        set_list_segments, 
        del_list_segments, 
        "list_segments's docstring")
    
    def insert_update_db(self, cur):
        for segment in self.__list_segments:
            segment.insert_update_db(cur)