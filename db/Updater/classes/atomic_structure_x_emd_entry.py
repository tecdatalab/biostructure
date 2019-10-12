'''
Created on 22 feb. 2019

@author: luis98
'''
from psycopg2 import sql


class Atomic_structure_x_emd_entry(object):
    '''
    classdocs
    '''
    __atomic_structure_entry_id = None
    __emd_entry_id = None

    def __init__(self, atomic_structure_entry_id, emd_entry_id):
        self.__atomic_structure_entry_id = atomic_structure_entry_id
        self.__emd_entry_id = emd_entry_id

    def insert_db(self, cur):
        cur.execute(
            sql.SQL("INSERT INTO atomic_structure_x_emd_entry(atomic_structure_entry_id,emd_entry_id) VALUES (%s,%s);"), [
                self.__atomic_structure_entry_id, self.__emd_entry_id])

    def get_atomic_structure_entry_id(self):
        return self.__atomic_structure_entry_id

    def get_emd_entry_id(self):
        return self.__emd_entry_id

    def set_atomic_structure_entry_id(self, value):
        self.__atomic_structure_entry_id = value

    def set_emd_entry_id(self, value):
        self.__emd_entry_id = value

    def del_atomic_structure_entry_id(self):
        del self.__atomic_structure_entry_id

    def del_emd_entry_id(self):
        del self.__emd_entry_id

    def __eq__(self, atomic_structure_x_emd_entry):
        return self.atomic_structure_entry_id == atomic_structure_x_emd_entry.atomic_structure_entry_id \
            and self.emd_entry_id == atomic_structure_x_emd_entry.emd_entry_id

    atomic_structure_entry_id = property(
        get_atomic_structure_entry_id,
        set_atomic_structure_entry_id,
        del_atomic_structure_entry_id,
        "atomic_structure_entry_id's docstring")
    emd_entry_id = property(
        get_emd_entry_id,
        set_emd_entry_id,
        del_emd_entry_id,
        "emd_entry_id's docstring")
