'''
Created on 22 feb. 2019

@author: luis98
'''

class Emd_entry(object):
    '''
    classdocs
    '''
    __id = None
    __full_name = None
    __acronym = None
    __volume = None
    __resolution = None
    __image_url = None
    __xml_url = None
    __map_url = None
    __map_information_id = None

    def __init__(self, id, full_name, acronym, volume, resolution, image_url, xml_url, map_url, map_information_id):
        self.__id = id
        self.__full_name = full_name
        self.__acronym = acronym
        self.__volume = volume
        self.__resolution = resolution
        self.__image_url = image_url
        self.__xml_url = xml_url
        self.__map_url = map_url
        self.__map_information_id = map_information_id

    def get_id(self):
        return self.__id


    def get_full_name(self):
        return self.__full_name


    def get_acronym(self):
        return self.__acronym


    def get_volume(self):
        return self.__volume


    def get_resolution(self):
        return self.__resolution


    def get_image_url(self):
        return self.__image_url


    def get_xml_url(self):
        return self.__xml_url


    def get_map_url(self):
        return self.__map_url


    def get_map_information_id(self):
        return self.__map_information_id


    def set_id(self, value):
        self.__id = value


    def set_full_name(self, value):
        self.__full_name = value


    def set_acronym(self, value):
        self.__acronym = value


    def set_volume(self, value):
        self.__volume = value


    def set_resolution(self, value):
        self.__resolution = value


    def set_image_url(self, value):
        self.__image_url = value


    def set_xml_url(self, value):
        self.__xml_url = value


    def set_map_url(self, value):
        self.__map_url = value


    def set_map_information_id(self, value):
        self.__map_information_id = value


    def del_id(self):
        del self.__id


    def del_full_name(self):
        del self.__full_name


    def del_acronym(self):
        del self.__acronym


    def del_volume(self):
        del self.__volume


    def del_resolution(self):
        del self.__resolution


    def del_image_url(self):
        del self.__image_url


    def del_xml_url(self):
        del self.__xml_url


    def del_map_url(self):
        del self.__map_url


    def del_map_information_id(self):
        del self.__map_information_id

    id = property(get_id, set_id, del_id, "id's docstring")
    full_name = property(get_full_name, set_full_name, del_full_name, "full_name's docstring")
    acronym = property(get_acronym, set_acronym, del_acronym, "acronym's docstring")
    volume = property(get_volume, set_volume, del_volume, "volume's docstring")
    resolution = property(get_resolution, set_resolution, del_resolution, "resolution's docstring")
    image_url = property(get_image_url, set_image_url, del_image_url, "image_url's docstring")
    xml_url = property(get_xml_url, set_xml_url, del_xml_url, "xml_url's docstring")
    map_url = property(get_map_url, set_map_url, del_map_url, "map_url's docstring")
    map_information_id = property(get_map_information_id, set_map_information_id, del_map_information_id, "map_information_id's docstring")


