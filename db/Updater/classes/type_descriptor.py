'''
Created on 22 feb. 2019

@author: luis98
'''
from psycopg2 import sql


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

    def insert_db(self, cur):
        if self.__id is None:
            cur.execute(
                sql.SQL("INSERT INTO type_descriptor(id,name,description) VALUES (DEFAULT,%s,%s) RETURNING id;"), [
                    self.__name, self.__description])
            self.id = [record for record in cur][0]
        else:
            cur.execute(
                sql.SQL("INSERT INTO type_descriptor(id,name,description) VALUES (%s,%s,%s);"), [
                    self.__id, self.__name, self.__description])

    def update_db(self, cur):
        cur.execute(
            sql.SQL("UPDATE type_descriptor SET name = %s , description = %s WHERE id = %s;"), [
                self.__name, self.__description, self.__id])

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

    def __eq__(self, type_descriptor):
        return self.id == type_descriptor.id \
            and self.name == type_descriptor.name \
            and self.description == type_descriptor.description

    id = property(get_id, set_id, del_id, "id's docstring")
    name = property(get_name, set_name, del_name, "name's docstring")
    description = property(
        get_description,
        set_description,
        del_description,
        "description's docstring")
