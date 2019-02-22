'''
Created on 22 feb. 2019

@author: luis98
'''

class Type_descriptor(object):
    '''
    classdocs
    '''
    __id = None
    __name = None
    __description = None

    def __init__(self, id, name, description):
        self.__id = id
        self.__name = name
        self.__description = description

    def get_id(self):
        return self.__id


    def get_name(self):
        return self.__name


    def get_description(self):
        return self.__description


    def set_id(self, value):
        self.__id = value


    def set_name(self, value):
        self.__name = value


    def set_description(self, value):
        self.__description = value


    def del_id(self):
        del self.__id


    def del_name(self):
        del self.__name


    def del_description(self):
        del self.__description

    id = property(get_id, set_id, del_id, "id's docstring")
    name = property(get_name, set_name, del_name, "name's docstring")
    description = property(get_description, set_description, del_description, "description's docstring")

    
    