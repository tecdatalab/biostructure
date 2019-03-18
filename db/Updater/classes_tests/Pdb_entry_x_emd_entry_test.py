'''
Created on 18 mar. 2019

@author: luis98
'''

from classes.Pdb_entry_x_emd_entry import Pdb_entry_x_emd_entry
import pytest
from pytest_postgresql import factories # importa la habilidad de usar postgresql

pdb_entry_x_emd_entry = Pdb_entry_x_emd_entry(None,None)

def test_0_set_get_id():
    
    for i in range(4):    
        pdb_entry_id = i            
        pdb_entry_x_emd_entry.set_pdb_entry_id(pdb_entry_id)
        assert pdb_entry_x_emd_entry.get_pdb_entry_id() == pdb_entry_id

def test_1_set_get_id():
    
    for i in range(4):    
        emd_entry_id = i            
        pdb_entry_x_emd_entry.set_emd_entry_id(emd_entry_id)
        assert pdb_entry_x_emd_entry.get_emd_entry_id() == emd_entry_id
        
def test_3_insert_db(postgresql):
    
    create_pdb_entry_x_emd_entry_table(postgresql)
    pdb_entry_x_emd_entrys = create_pdb_entry_x_emd_entrys()
    insert_pdb_entry_x_emd_entrys(postgresql,pdb_entry_x_emd_entrys)
    result = get_entry_x_emd_entrys(postgresql)
    assert result[0] == Pdb_entry_x_emd_entry(1,1)
  
def get_entry_x_emd_entrys(connection):
    
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM pdb_entry_x_emd_entry")
    type_descriptors = [Pdb_entry_x_emd_entry(record[0],record[1]) for record in cursor]
    cursor.close()

    return type_descriptors  
    
def insert_pdb_entry_x_emd_entrys(connection, pdb_entry_x_emd_entrys):
    cursor = connection.cursor()
    
    for pdb_entry_x_emd_entry in pdb_entry_x_emd_entrys:
        pdb_entry_x_emd_entry.insert_db(cursor)

    connection.commit()
    cursor.close()      
    
def create_pdb_entry_x_emd_entrys():
    resultado = []
    resultado.append(Pdb_entry_x_emd_entry(1,1))
    return resultado
    
def create_pdb_entry_x_emd_entry_table(connection):
    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE pdb_entry_x_emd_entry(\
        pdb_entry_id INT NOT NULL, \
        emd_entry_id INT NOT NULL, \
        PRIMARY KEY (emd_entry_id,pdb_entry_id)\
        );")
    connection.commit()
    cursor.close()
        