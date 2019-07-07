'''
Created on 22 feb. 2019

@author: luis98
'''
import time
import datetime
from datetime import date
from psycopg2 import sql

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
        if isinstance(modification, datetime.date):
            self.__modification = modification
        elif modification == None:
            self.__modification = None
        else:
            self.__modification = date(modification[0],modification[1],modification[2])
        self.__map_file = map_file
        self.__xml_file = xml_file
        self.__image_file = image_file
        
    def insert_db(self, cur):
        cur.execute(sql.SQL("INSERT INTO time_stamp(emd_entry_id,modification,map_file,xml_file,image_file) VALUES (%s,%s,%s,%s,%s);")
        ,[self.__emd_entry_id,self.__modification,self.__map_file,self.__xml_file,self.__image_file])
        
        
    def update_db(self, cur):
        cur.execute(sql.SQL("UPDATE time_stamp set modification = %s, map_file = %s, xml_file = %s, image_file = %s WHERE emd_entry_id = %s;")
        ,[self.__modification,self.__map_file,self.__xml_file,self.__image_file,self.__emd_entry_id])


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
        if isinstance(value, datetime.date):
            self.__modification = value
        elif value == None:
            self.__modification = None
        else:
            self.__modification = date(value[0],value[1],value[2])

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
        
        
    def __eq__(self, time_stamp):        
        return self.emd_entry_id == time_stamp.emd_entry_id \
           and self.modification == time_stamp.modification \
           and self.map_file == time_stamp.map_file \
           and self.xml_file == time_stamp.xml_file \
           and self.image_file == time_stamp.image_file

    emd_entry_id = property(get_emd_entry_id, set_emd_entry_id, del_emd_entry_id, "emd_entry_id's docstring")
    modification = property(get_modification, set_modification, del_modification, "modification's docstring")
    map_file = property(get_map_file, set_map_file, del_map_file, "map_file's docstring")
    xml_file = property(get_xml_file, set_xml_file, del_xml_file, "xml_file's docstring")
    image_file = property(get_image_file, set_image_file, del_image_file, "image_file's docstring")
        
    

