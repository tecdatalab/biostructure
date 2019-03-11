'''
Created on 3 mar. 2019

@author: luis98
'''
from classes.Type_descriptor import Type_descriptor
import pytest
from pytest_postgresql import factories # importa la habilidad de usar postgresql

type_descriptor = Type_descriptor(None,None,None)


def test_0_set_get_id():
    
    for i in range(4):    
        id = i            
        type_descriptor.set_id(id)
        assert type_descriptor.get_id() == id

def test_1_set_get_name():
    
    for i in range(4):    
        name = 'name' + str(i)            
        type_descriptor.set_name(name)
        assert type_descriptor.get_name() == name

def test_2_set_get_description():
    
    for i in range(4):    
        description = 'description' + str(i)            
        type_descriptor.set_description(description)
        assert type_descriptor.get_description() == description
        
def test_3_insert_db(postgresql):
    
    create_type_descriptor_table(postgresql)
    types_descriptors = create_type_descriptors()
    insert_type_descriptors(postgresql,types_descriptors)
    result = get_type_descriptors(postgresql)
    for i in range(len(result)):
        assert result[i] == types_descriptors[i]
        
def test_4_update_db(postgresql):
    
    create_type_descriptor_table(postgresql)
    types_descriptors = create_type_descriptors()
    insert_type_descriptors(postgresql,types_descriptors)
    
    for i in types_descriptors:
        i.name = i.name + " update"
        i.description = i.description + " update"
    
    update_type_descriptors(postgresql, types_descriptors)
    result = get_type_descriptors(postgresql)
    for i in range(len(result)):
        assert result[i] == types_descriptors[i]
    
        
def create_type_descriptor_table(connection):
    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE type_descriptor(\
        id SERIAL PRIMARY KEY,\
        name TEXT UNIQUE NOT NULL,\
        description TEXT NOT NULL\
        );")
    connection.commit()
    cursor.close()
    
def create_type_descriptors():
    resultado = []
    for i in range(5):
        resultado.append(Type_descriptor(i+1,"Name "+str(i) ,"Descripcion "+str(i)))
    return resultado

def insert_type_descriptors(connection, type_descriptors):
    cursor = connection.cursor()
    
    for type_descriptor in type_descriptors:
        type_descriptor.insert_db(cursor)

    connection.commit()
    cursor.close()
    
def update_type_descriptors(connection, type_descriptors):
    cursor = connection.cursor()
    
    for type_descriptor in type_descriptors:
        type_descriptor.update_db(cursor)

    connection.commit()
    cursor.close()

def get_type_descriptors(connection):
    
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM type_descriptor")
    type_descriptors = [Type_descriptor(record[0],record[1],record[2]) for record in cursor]
    cursor.close()

    return type_descriptors


