'''
Created on 10 mar. 2019

@author: luis98
'''
from classes.Volume_filter import Volume_filter
import pytest
from pytest_postgresql import factories # Import the ability to use postgresql

volume_filter = Volume_filter(None,None)


def test_0_set_get_id():
    
    for i in range(4):    
        id = i            
        volume_filter.set_id(id)
        assert volume_filter.get_id() == id

def test_1_set_get_name():
    
    for i in range(4):    
        name = 'name' + str(i)            
        volume_filter.set_name(name)
        assert volume_filter.get_name() == name

def test_3_insert_db(postgresql):
    
    create_volume_filter_table(postgresql)
    volume_filters = create_volume_filters()
    insert_volume_filters(postgresql,volume_filters)
    result = get_volume_filters(postgresql)
    assert result[0] == Volume_filter(1,"Name 1")
    
def test_4_update_db(postgresql):
    
    create_volume_filter_table(postgresql)
    volume_filters = create_volume_filters()
    insert_volume_filters(postgresql,volume_filters)
    volume_filters[0].name = "Name 2"
    
    update_volume_filters(postgresql,volume_filters)
    result = get_volume_filters(postgresql)
    assert result[0] == Volume_filter(1,"Name 2")
    
def get_volume_filters(connection):
    
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM volume_filter")
    type_descriptors = [Volume_filter(record[0],record[1]) for record in cursor]
    cursor.close()

    return type_descriptors

def update_volume_filters(connection, volume_filters):
    cursor = connection.cursor()
    
    for volume_filter in volume_filters:
        volume_filter.update_db(cursor)

    connection.commit()
    cursor.close()  
    
def insert_volume_filters(connection, volume_filters):
    cursor = connection.cursor()
    
    for volume_filter in volume_filters:
        volume_filter.insert_db(cursor)

    connection.commit()
    cursor.close()    
    
def create_volume_filters():
    resultado = []
    resultado.append(Volume_filter(1,"Name 1"))
    return resultado

def create_volume_filter_table(connection):
    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE volume_filter(\
        id SERIAL PRIMARY KEY,\
        name TEXT UNIQUE NOT NULL\
        );")
    connection.commit()
    cursor.close()








