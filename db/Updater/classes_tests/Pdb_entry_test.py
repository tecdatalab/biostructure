'''
Created on 18 mar. 2019

@author: luis98
'''
from classes.Pdb_entry import Pdb_entry
import pytest
from pytest_postgresql import factories # Import the ability to use postgresql

pdb_entry = Pdb_entry(None,None)

def test_0_set_get_id():
    
    for i in range(4):    
        id = i            
        pdb_entry.set_id(id)
        assert pdb_entry.get_id() == id

def test_1_set_get_pdb():
    
    for i in range(4):    
        pdb = 'pdb' + str(i)            
        pdb_entry.set_pdb(pdb)
        assert pdb_entry.get_pdb() == pdb
        
def test_3_insert_db(postgresql):
    
    create_pdb_entry_table(postgresql)
    pdb_entrys = create_pdb_entrys()
    insert_pdb_entrys(postgresql,pdb_entrys)
    result = get_pdb_entrys(postgresql)
    assert result[0] == Pdb_entry(1,"pbd1")
    
def test_4_update_db(postgresql):
    
    create_pdb_entry_table(postgresql)
    pdb_entrys = create_pdb_entrys()
    insert_pdb_entrys(postgresql,pdb_entrys)
    pdb_entrys[0].pdb = "pbd1Up"
    update_pdb_entrys(postgresql,pdb_entrys)
    
    result = get_pdb_entrys(postgresql)
    assert result[0] == Pdb_entry(1,"pbd1Up")

def get_pdb_entrys(connection):
    
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM pdb_entry")
    type_descriptors = [Pdb_entry(record[0],record[1]) for record in cursor]
    cursor.close()

    return type_descriptors

def update_pdb_entrys(connection, pdb_entrys):
    cursor = connection.cursor()
    
    for pdb_entry in pdb_entrys:
        pdb_entry.update_db(cursor)

    connection.commit()
    cursor.close()  

def insert_pdb_entrys(connection, pdb_entrys):
    cursor = connection.cursor()
    
    for pdb_entry in pdb_entrys:
        pdb_entry.insert_db(cursor)

    connection.commit()
    cursor.close()  

def create_pdb_entrys():
    resultado = []
    resultado.append(Pdb_entry(1,"pbd1"))
    return resultado
        
def create_pdb_entry_table(connection):
    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE pdb_entry(\
        id SERIAL PRIMARY KEY,\
        pdb TEXT NOT NULL\
        );")
    connection.commit()
    cursor.close()
        