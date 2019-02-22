'''
Created on 22 feb. 2019

@author: luis98
'''

class Time_stamp(object):
    '''
    classdocs
    '''
    __emd_entry_id = None
    __modification = None
    __map_file = None
    __xml_file = None
    __image_file = None

    def __init__(self, emd_entry_id, modification, map_file, xml_file, image_file):
        self.__emd_entry_id = emd_entry_id
        self.__modification = modification
        self.__map_file = map_file
        self.__xml_file = xml_file
        self.__image_file = image_file

    def get_emd_entry_id(self):
        return self.__emd_entry_id


    def get_modification(self):
        return self.__modification


    def get_map_file(self):
        return self.__map_file


    def get_xml_file(self):
        return self.__xml_file


    def get_image_file(self):
        return self.__image_file


    def set_emd_entry_id(self, value):
        self.__emd_entry_id = value


    def set_modification(self, value):
        self.__modification = value


    def set_map_file(self, value):
        self.__map_file = value


    def set_xml_file(self, value):
        self.__xml_file = value


    def set_image_file(self, value):
        self.__image_file = value


    def del_emd_entry_id(self):
        del self.__emd_entry_id


    def del_modification(self):
        del self.__modification


    def del_map_file(self):
        del self.__map_file


    def del_xml_file(self):
        del self.__xml_file


    def del_image_file(self):
        del self.__image_file

    emd_entry_id = property(get_emd_entry_id, set_emd_entry_id, del_emd_entry_id, "emd_entry_id's docstring")
    modification = property(get_modification, set_modification, del_modification, "modification's docstring")
    map_file = property(get_map_file, set_map_file, del_map_file, "map_file's docstring")
    xml_file = property(get_xml_file, set_xml_file, del_xml_file, "xml_file's docstring")
    image_file = property(get_image_file, set_image_file, del_image_file, "image_file's docstring")
        
    

