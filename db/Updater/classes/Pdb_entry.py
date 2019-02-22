'''
Created on 22 feb. 2019

@author: luis98
'''

class Pdb_entry(object):
    '''
    classdocs
    '''
    __id = None
    __pdb = None

    def __init__(self, id, pdb):
        self.__id = id
        self.__pdb = pdb

    def get_id(self):
        return self.__id


    def get_pdb(self):
        return self.__pdb


    def set_id(self, value):
        self.__id = value


    def set_pdb(self, value):
        self.__pdb = value


    def del_id(self):
        del self.__id


    def del_pdb(self):
        del self.__pdb

    id = property(get_id, set_id, del_id, "id's docstring")
    pdb = property(get_pdb, set_pdb, del_pdb, "pdb's docstring")

        