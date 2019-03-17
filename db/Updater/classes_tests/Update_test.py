'''
Created on 17 mar. 2019

@author: luis98
'''
from classes.Update import Update
import pytest
from pytest_postgresql import factories # importa la habilidad de usar postgresql
import time
from datetime import date

update = Update(None)

def test_0_set_get_last_update():
       
    update.set_last_update((2018,2,14))
    assert update.get_last_update() == date(2018,2,14)

def test_1_insert_db(postgresql):
    
    create_update_table(postgresql)
    updates = create_updates()
    insert_updates(postgresql,updates)
    result = get_volume_filters(postgresql)
    assert result[0] == Update((2018,2,14))

def get_volume_filters(connection):
    
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM update")
    type_descriptors = [Update(record[0]) for record in cursor]
    cursor.close()

    return type_descriptors

def insert_updates(connection, updates):
    cursor = connection.cursor()
    
    for update in updates:
        update.insert_db(cursor)

    connection.commit()
    cursor.close()    

def create_updates():
    resultado = []
    resultado.append(Update((2018,2,14)))
    return resultado

def create_update_table(connection):
    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE update(\
        last_update DATE PRIMARY KEY\
        );")
    connection.commit()
    cursor.close()

