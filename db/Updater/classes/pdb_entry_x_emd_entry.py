'''
Created on 22 feb. 2019

@author: luis98
'''
from psycopg2 import sql


class Pdb_entry_x_emd_entry(object):
    '''
    classdocs
    '''
    __pdb_entry_id = None
    __emd_entry_id = None

    def __init__(self, pdb_entry_id, emd_entry_id):
        self.__pdb_entry_id = pdb_entry_id
        self.__emd_entry_id = emd_entry_id

    def insert_db(self, cur):
        cur.execute(
            sql.SQL("INSERT INTO pdb_entry_x_emd_entry(pdb_entry_id,emd_entry_id) VALUES (%s,%s);"), [
                self.__pdb_entry_id, self.__emd_entry_id])
    
    def insert_update_db(self, cur):
        cur.execute(
            sql.SQL("INSERT INTO pdb_entry_x_emd_entry(pdb_entry_id,emd_entry_id) VALUES (%s,%s);"), [
                self.__pdb_entry_id, self.__emd_entry_id])

    def get_pdb_entry_id(self):
        return self.__pdb_entry_id

    def get_emd_entry_id(self):
        return self.__emd_entry_id

    def set_pdb_entry_id(self, value):
        self.__pdb_entry_id = value

    def set_emd_entry_id(self, value):
        self.__emd_entry_id = value

    def del_pdb_entry_id(self):
        del self.__pdb_entry_id

    def del_emd_entry_id(self):
        del self.__emd_entry_id

    def __eq__(self, pdb_entry_x_emd_entry):
        return self.pdb_entry_id == pdb_entry_x_emd_entry.pdb_entry_id \
            and self.emd_entry_id == pdb_entry_x_emd_entry.emd_entry_id

    pdb_entry_id = property(
        get_pdb_entry_id,
        set_pdb_entry_id,
        del_pdb_entry_id,
        "pdb_entry_id's docstring")
    emd_entry_id = property(
        get_emd_entry_id,
        set_emd_entry_id,
        del_emd_entry_id,
        "emd_entry_id's docstring")
