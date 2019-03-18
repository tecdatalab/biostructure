'''
Created on 10 mar. 2019

@author: luis98
'''

from classes.Descriptor import Descriptor
import pytest
from pytest_postgresql import factories # Import the ability to use postgresql

descriptor = Descriptor(None,None,None)

def test_0_set_get_emd_entry_id():
    
    for i in range(4):    
        emd_entry_id = i            
        descriptor.emd_entry_id = emd_entry_id
        assert descriptor.get_emd_entry_id() == emd_entry_id

def test_1_set_get_type_descriptor_id():
    
    for i in range(4):    
        type_descriptor_id = i            
        descriptor.type_descriptor_id = type_descriptor_id
        assert descriptor.get_type_descriptor_id() == type_descriptor_id

def test_2_set_get_numbers():
    
    for i in range(4):    
        numbers = list(range(-2, i))           
        descriptor.set_numbers(numbers)
        assert descriptor.get_numbers() == numbers
        
def test_3_insert_db(postgresql):
    
    create_descriptor_table(postgresql)
    descriptors = create_descriptors()
    insert_descriptors(postgresql,descriptors)
    result = get_descriptors(postgresql)
    assert result[0] == Descriptor(1,1,[1,2,3])

def test_4_update_db(postgresql):
    
    create_descriptor_table(postgresql)
    descriptors = create_descriptors()
    insert_descriptors(postgresql,descriptors)
    
    descriptors[0].numbers = [1,2,8,3,2]
    
    update_descriptors(postgresql, descriptors)
    result = get_descriptors(postgresql)
    assert result[0] == Descriptor(1,1,[1,2,8,3,2])

def create_descriptor_table(connection):
    
    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE descriptor(\
        emd_entry_id INT NOT NULL,\
        type_descriptor_id INT NOT NULL,\
        numbers JSON NOT NULL,\
        PRIMARY KEY (emd_entry_id,type_descriptor_id)\
        );")
    connection.commit()
    cursor.close()

def create_descriptors():
    resultado = []
    resultado.append(Descriptor(1,1,[1,2,3]))
    return resultado

def insert_descriptors(connection, descriptors):
    cursor = connection.cursor()
    
    for descriptor in descriptors:
        descriptor.insert_db(cursor)

    connection.commit()
    cursor.close()

def update_descriptors(connection, descriptors):
    cursor = connection.cursor()
    
    for descriptor in descriptors:
        descriptor.update_db(cursor)

    connection.commit()
    cursor.close()
    
def get_descriptors(connection):
    
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM descriptor")
    descriptors = [Descriptor(record[0],record[1],record[2]) for record in cursor]
    cursor.close()

    return descriptors
    
    


