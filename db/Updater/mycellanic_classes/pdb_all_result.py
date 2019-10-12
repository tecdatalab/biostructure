'''
Created on 22 sep. 2019

@author: luis98
'''


class Pdb_all_result:
    '''
    classdocs
    '''
    __id_code = None
    __structure_type = None
    __method = None
    

    def __init__(self, id_code, method, structure_type):
        self.__id_code = id_code
        self.__method = method
        self.__structure_type = structure_type

    def get_id_code(self):
        return self.__id_code

    def get_structure_type(self):
        return self.__structure_type
    
    def get_method(self):
        return self.__method

    def set_id_code(self, value):
        self.__id_code = value

    def set_structure_type(self, value):
        self.__structure_type = value
    
    def set_method(self, value):
        self.__method = value

    def del_id_code(self):
        del self.__id_code

    def del_structure_type(self):
        del self.__structure_type

    def del_method(self):
        del self.__method

    id_code = property(get_id_code, set_id_code, del_id_code, "id_code's docstring")
    structure_type = property(get_structure_type, set_structure_type, del_structure_type, "structure_type's docstring")
    method = property(get_method, set_method, del_method, "method's docstring")
    
