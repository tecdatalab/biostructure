'''
Created on 22 feb. 2019

@author: luis98
'''

class Search_history(object):
    '''
    classdocs
    '''
    __id = None
    __date_time = None
    __ip = None
    __emd_entry_id = None
    __name_file = None
    __counter_level = None
    __representation_id = None
    __volume_filter_id = None
    __resolution_filter_min = None
    __resolution_filter_max = None
    
    def __init__(self, id, date_time, ip, emd_entry_id, name_file, counter_level, representation_id, volume_filter_id, resolution_filter_min, resolution_filter_max):
        self.__id = id
        self.__date_time = date_time
        self.__ip = ip
        self.__emd_entry_id = emd_entry_id
        self.__name_file = name_file
        self.__counter_level = counter_level
        self.__representation_id = representation_id
        self.__volume_filter_id = volume_filter_id
        self.__resolution_filter_min = resolution_filter_min
        self.__resolution_filter_max = resolution_filter_max

    def get_id(self):
        return self.__id


    def get_date_time(self):
        return self.__date_time


    def get_ip(self):
        return self.__ip


    def get_emd_entry_id(self):
        return self.__emd_entry_id


    def get_name_file(self):
        return self.__name_file


    def get_counter_level(self):
        return self.__counter_level


    def get_representation_id(self):
        return self.__representation_id


    def get_volume_filter_id(self):
        return self.__volume_filter_id


    def get_resolution_filter_min(self):
        return self.__resolution_filter_min


    def get_resolution_filter_max(self):
        return self.__resolution_filter_max


    def set_id(self, value):
        self.__id = value


    def set_date_time(self, value):
        self.__date_time = value


    def set_ip(self, value):
        self.__ip = value


    def set_emd_entry_id(self, value):
        self.__emd_entry_id = value


    def set_name_file(self, value):
        self.__name_file = value


    def set_counter_level(self, value):
        self.__counter_level = value


    def set_representation_id(self, value):
        self.__representation_id = value


    def set_volume_filter_id(self, value):
        self.__volume_filter_id = value


    def set_resolution_filter_min(self, value):
        self.__resolution_filter_min = value


    def set_resolution_filter_max(self, value):
        self.__resolution_filter_max = value


    def del_id(self):
        del self.__id


    def del_date_time(self):
        del self.__date_time


    def del_ip(self):
        del self.__ip


    def del_emd_entry_id(self):
        del self.__emd_entry_id


    def del_name_file(self):
        del self.__name_file


    def del_counter_level(self):
        del self.__counter_level


    def del_representation_id(self):
        del self.__representation_id


    def del_volume_filter_id(self):
        del self.__volume_filter_id


    def del_resolution_filter_min(self):
        del self.__resolution_filter_min


    def del_resolution_filter_max(self):
        del self.__resolution_filter_max

    id = property(get_id, set_id, del_id, "id's docstring")
    date_time = property(get_date_time, set_date_time, del_date_time, "date_time's docstring")
    ip = property(get_ip, set_ip, del_ip, "ip's docstring")
    emd_entry_id = property(get_emd_entry_id, set_emd_entry_id, del_emd_entry_id, "emd_entry_id's docstring")
    name_file = property(get_name_file, set_name_file, del_name_file, "name_file's docstring")
    counter_level = property(get_counter_level, set_counter_level, del_counter_level, "counter_level's docstring")
    representation_id = property(get_representation_id, set_representation_id, del_representation_id, "representation_id's docstring")
    volume_filter_id = property(get_volume_filter_id, set_volume_filter_id, del_volume_filter_id, "volume_filter_id's docstring")
    resolution_filter_min = property(get_resolution_filter_min, set_resolution_filter_min, del_resolution_filter_min, "resolution_filter_min's docstring")
    resolution_filter_max = property(get_resolution_filter_max, set_resolution_filter_max, del_resolution_filter_max, "resolution_filter_max's docstring")
        
    
