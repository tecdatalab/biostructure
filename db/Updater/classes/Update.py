'''
Created on 22 feb. 2019

@author: luis98
'''
import time
import datetime
from datetime import date
from psycopg2 import sql

'''Format YYYY-MM-DD'''
class Update(object):
    '''
    classdocs
    '''
    __last_update = None

    def __init__(self, last_update):
        if isinstance(last_update, datetime.date):
            self.__last_update = last_update
        elif last_update == None:
            self.__last_update = None
        else:
            self.__last_update = date(last_update[0],last_update[1],last_update[2])
            
    def insert_db(self, cur):
        cur.execute(sql.SQL("INSERT INTO update(last_update) VALUES (%s);")
        ,[self.__last_update])
        

    def get_last_update(self):
        return self.__last_update


    def set_last_update(self, value):
        if isinstance(value, datetime.date):
            self.__last_update = value
        elif value == None:
            self.__last_update = None
        else:
            self.__last_update = date(value[0],value[1],value[2])


    def del_last_update(self):
        del self.__last_update
    
    def __eq__(self, update):        
        return self.last_update == update.last_update
           
    last_update = property(get_last_update, set_last_update, del_last_update, "last_update's docstring")      