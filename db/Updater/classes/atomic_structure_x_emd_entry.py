'''
Created on 22 feb. 2019

@author: luis98
'''
from psycopg2 import sql


class Atomic_structure_x_emd_entry(object):
    '''
    classdocs
    '''
    __atomic_structure = None
    __emd_entry = None
    __atomic_structure_id = None
    __emd_entry_id = None


    def __init__(self, atomic_structure, emd_entry):
        self.__atomic_structure = atomic_structure
        self.__emd_entry = int(emd_entry)
        self.__atomic_structure_id = None
        
    def __insert_db(self, cur):
        cur.execute(
            sql.SQL("INSERT INTO atomic_structure_x_emd_entry(atomic_structure,emd_entry) VALUES (%s,%s);"), [
                self.__atomic_structure_id, self.__emd_entry])

    def insert_update_db(self, cur):
        if self.__atomic_structure_id == None:
            cur.execute(
            sql.SQL("SELECT id FROM atomic_structure WHERE id_code = %s;"), [
                self.__atomic_structure])
            result = [record[0] for record in cur]
            self.__atomic_structure_id = result[0]           

        cur.execute(
            sql.SQL("SELECT * FROM atomic_structure_x_emd_entry WHERE atomic_structure = %s AND emd_entry = %s;"), [
                self.__atomic_structure_id, self.__emd_entry])
        result = [record[0] for record in cur]
        if len(result)==0:
            self.__insert_db(cur)
            
    def get_atomic_structure(self):
        return self.__atomic_structure

    def get_emd_entry(self):
        return self.__emd_entry

    def set_atomic_structure(self, value):
        self.__atomic_structure = value

    def set_emd_entry(self, value):
        self.__emd_entry = int(value)

    def del_atomic_structure(self):
        del self.__atomic_structure

    def del_emd_entry(self):
        del self.__emd_entry

    def __eq__(self, atomic_structure_x_emd_entry):
        return self.atomic_structure == atomic_structure_x_emd_entry.atomic_structure \
            and self.emd_entry == atomic_structure_x_emd_entry.emd_entry

    atomic_structure = property(
        get_atomic_structure,
        set_atomic_structure,
        del_atomic_structure,
        "atomic_structure's docstring")
    emd_entry = property(
        get_emd_entry,
        set_emd_entry,
        del_emd_entry,
        "emd_entry's docstring")
