'''
Created on 29 jun. 2020

@author: dnnxl
'''
import datetime
from datetime import date
from psycopg2 import sql

'''Format YYYY-MM-DD'''

class Binnacle(object):
    '''
    classdocs
    '''
    __last_update = None
    __emd_id = None
    __pdb_id = None

    def __init__(self):
        self.__last_update = None

    def set_emd_id(self, emd_id):
        self.__emd_id = emd_id
    
    def set_pdb_id(self, pdb_id):
        self.__pdb_id = pdb_id
    
    def get_pdb_id(self):
        return self.__pdb_id

    def get_emd_id(self):
        return self.__emd_id
    
    def get_last_update(self):
        return self.__last_update
    
    def set_last_update(self, cur):
        cur.execute("SELECT * FROM updater ORDER BY last_update DESC")
        updates = [record[0] for record in cur]
        self.__last_update = updates[0]

    def update_binnacle_emd(self, cur):
        cur.execute(
            sql.SQL('CALL update_binnacle_emd({0}, \'{1}\', {2});'.format(self.__emd_id, self.__last_update, 1)))

    def update_binnacle_pdb(self, cur):
        cur.execute(
            sql.SQL('CALL update_binnacle_pdb(\'{0}\', \'{1}\', {2});'.format(self.__pdb_id, self.__last_update, 1)))

    def insert_binnacle_emd(self, cur):
        cur.execute(
            sql.SQL('CALL insert_binnacle_emd({0}, \'{1}\', {2});'.format(self.__emd_id, self.__last_update, 1)))

    def insert_binnacle_pdb(self, cur):
        cur.execute(
            sql.SQL('CALL insert_binnacle_pdb(\'{0}\', \'{1}\', {2});'.format(self.__pdb_id, self.__last_update, 1)))

    def get_attempt_pdb(self, cur):
        cur.execute(
            sql.SQL('SELECT get_attempt_pdb(\'{0}\', \'{1}\');'.format(self.__pdb_id, self.__last_update)))
        return cur.fetchone()[0]

    def get_attempt_emd(self, cur):
        cur.execute(
            sql.SQL('SELECT get_attempt_emd({0}, \'{1}\');'.format(self.__emd_id, self.__last_update)))
        return cur.fetchone()[0]

    