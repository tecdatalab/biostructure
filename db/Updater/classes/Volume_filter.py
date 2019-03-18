'''
Created on 22 feb. 2019

@author: luis98
'''
from psycopg2 import sql

class Volume_filter(object):
    '''
    classdocs
    '''
    __id = None
    __name = None
    
    def __init__(self, id, name):
        self.__id = id
        self.__name = name
        
        
    def insert_db(self, cur):
        cur.execute(sql.SQL("INSERT INTO volume_filter(id,name) VALUES (DEFAULT,%s);")
        ,[self.__name])
        

    def update_db(self, cur):
        cur.execute(sql.SQL("UPDATE volume_filter SET name = %s WHERE id = %s;")
        ,[self.__name,self.__id])


    def get_id(self):
        return self.__id


    def get_name(self):
        return self.__name


    def set_id(self, value):
        self.__id = value


    def set_name(self, value):
        self.__name = value


    def del_id(self):
        del self.__id


    def del_name(self):
        del self.__name
    
        
    def __eq__(self, volume_filter):        
        return self.id == volume_filter.id \
           and self.name == volume_filter.name

    id = property(get_id, set_id, del_id, "id's docstring")
    name = property(get_name, set_name, del_name, "name's docstring")

