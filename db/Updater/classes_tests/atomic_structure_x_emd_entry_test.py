'''
Created on 27 oct. 2019

@author: luis98
'''
from classes.atomic_structure_x_emd_entry import Atomic_structure_x_emd_entry 
from classes.atomic_structure import Atomic_structure
from classes.emd_entry import Emd_entry
import pytest
from pytest_postgresql import factories  # Import the ability to use postgresql

atomic_structure_x_emd_entry = Atomic_structure_x_emd_entry(None, 0)


def test_0_set_get_atomic_structure():
    atomic_structure_x_emd_entry.atomic_structure = '100d'
    assert atomic_structure_x_emd_entry.atomic_structure == '100d'

def test_1_set_get_emd_entry():
    atomic_structure_x_emd_entry.emd_entry = 1
    assert atomic_structure_x_emd_entry.emd_entry == 1

def test_3_insert_db(postgresql):
    create_atomic_structure_x_emd_entry_table(postgresql)
    atomic_structure_x_emd_entrys = create_atomic_structure_x_emd_entry(postgresql)
    insert_atomic_structure_x_emd_entry(postgresql, atomic_structure_x_emd_entrys)
    result = get_atomic_structure_x_emd_entry(postgresql)
    assert result[0].atomic_structure == 1 and result[0].emd_entry == 1

def test_4_insert_update_db(postgresql):
    create_atomic_structure_x_emd_entry_table(postgresql)
    atomic_structure_x_emd_entrys = create_atomic_structure_x_emd_entry(postgresql)
    insert_atomic_structure_x_emd_entry(postgresql, atomic_structure_x_emd_entrys)
    insert_update_atomic_structure_x_emd_entry(postgresql, atomic_structure_x_emd_entrys)
    result = get_atomic_structure_x_emd_entry(postgresql)
    assert result[0].atomic_structure == 1 and result[0].emd_entry == 1
 
def create_atomic_structure_x_emd_entry(connection):
    cursor = connection.cursor()
    temp_at = Atomic_structure(None, '100d', None, 1, None, None, None, None, None, None, None)
    temp_em = Emd_entry(id = 1)
    temp_at.insert_update_db_complex(cursor)
    temp_em.insert_update_db(cursor)
    result = []
    result.append(Atomic_structure_x_emd_entry('100d', 1))
    return result

def insert_atomic_structure_x_emd_entry(connection, atomic_structure_x_emd_entrys):
    cursor = connection.cursor()

    for atomic_structure_x_emd_entry in atomic_structure_x_emd_entrys:
        atomic_structure_x_emd_entry.insert_update_db(cursor)

    connection.commit()
    cursor.close()

def insert_update_atomic_structure_x_emd_entry(connection, atomic_structure_x_emd_entrys):
    cursor = connection.cursor()

    for atomic_structure_x_emd_entry in atomic_structure_x_emd_entrys:
        atomic_structure_x_emd_entry.insert_update_db(cursor)

    connection.commit()
    cursor.close()
    
def get_atomic_structure_x_emd_entry(connection):

    cursor = connection.cursor()
    cursor.execute("SELECT * FROM atomic_structure_x_emd_entry")
    atomic_structure_x_emd_entrys = [Atomic_structure_x_emd_entry(record[0], record[1])
                   for record in cursor]
    cursor.close()

    return atomic_structure_x_emd_entrys


def create_atomic_structure_x_emd_entry_table(connection):

    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE atomic_structure_x_emd_entry(\
        atomic_structure_id INT,\
        emd_entry_id INT\
        );")
    
    cursor.execute(
        "CREATE TABLE atomic_structure(\
        id SERIAL PRIMARY KEY,\
        id_code TEXT UNIQUE,\
        parent INT,\
        atomic_structure_type_id INT,\
        sequence TEXT,\
        atoms TEXT,\
        atoms_count INT,\
        aminoacid_count INT,\
        png_img_3d TEXT,\
        gif_img_3d TEXT,\
        numbers_descriptor JSON\
        );")
    
    cursor.execute(
        "CREATE TABLE emd_entry(\
        id INT PRIMARY KEY,\
        full_name TEXT,\
        acronym TEXT,\
        volume FLOAT8,\
        resolution FLOAT8,\
        image_url TEXT,\
        png_img_3d TEXT,\
        gif_img_3d TEXT,\
        xml_url TEXT,\
        map_url TEXT,\
        map_information_id INT\
        );")
    
    cursor.execute(
        "CREATE TABLE map_information(\
        id SERIAL PRIMARY KEY,\
        file_information JSON,\
        data_type TEXT,\
        num_columns INT,\
        num_rows INT,\
        num_sections INT,\
        origin_col INT,\
        origin_row INT,\
        origin_sec INT,\
        limit_col INT,\
        limit_row INT,\
        limit_sec INT,\
        spacing_col INT,\
        spacing_row INT,\
        spacing_sec INT,\
        cell_a JSON,\
        cell_b JSON,\
        cell_c JSON,\
        cell_alpha JSON,\
        cell_beta JSON,\
        cell_gamma JSON,\
        axis_order_fast CHAR(1),\
        axis_order_medium CHAR(1),\
        axis_order_slow CHAR(1),\
        minimum FLOAT8,\
        maximum FLOAT8,\
        average FLOAT8,\
        std FLOAT8,\
        space_group_number INT,\
        details TEXT,\
        pixel_x JSON,\
        pixel_y JSON,\
        pixel_z JSON,\
        countour_level FLOAT8,\
        annotation_details TEXT\
        );")
    
    connection.commit()
    cursor.close()