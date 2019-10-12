'''
Created on 17 mar. 2019

@author: luis98
'''
from classes.Time_stamp import Time_stamp
import pytest
from pytest_postgresql import factories  # Import the ability to use postgresql
import time
from datetime import date

update = Time_stamp(None, None, None, None, None)


def test_0_set_get_emd_entry_id():

    update.set_emd_entry_id(1)
    assert update.get_emd_entry_id() == 1


def test_1_set_get_modification():

    update.set_modification((2018, 2, 14))
    assert update.get_modification() == date(2018, 2, 14)


def test_2_set_get_map_file():

    update.set_map_file("mapFile")
    assert update.get_map_file() == "mapFile"


def test_3_set_get_xml_file():

    update.set_xml_file("xmlFile")
    assert update.get_xml_file() == "xmlFile"


def test_4_set_get_image_file():

    update.set_image_file("imageFile")
    assert update.get_image_file() == "imageFile"


def test_5_insert_db(postgresql):

    create_time_stamp_table(postgresql)
    time_stamps = create_time_stamps()
    insert_time_stamps(postgresql, time_stamps)
    result = get_time_stamps(postgresql)
    assert result[0] == Time_stamp(
        1, (2018, 2, 14), "mapFile", "xmlFile", "imageFile")


def test_6_update_db(postgresql):

    create_time_stamp_table(postgresql)
    time_stamps = create_time_stamps()
    insert_time_stamps(postgresql, time_stamps)

    time_stamps[0].modification = (2018, 2, 15)
    time_stamps[0].map_file = "mapFileup"
    time_stamps[0].xml_file = "xmlFileup"
    time_stamps[0].image_file = "imageFileup"

    update_time_stamps(postgresql, time_stamps)
    result = get_time_stamps(postgresql)
    assert result[0] == Time_stamp(
        1, (2018, 2, 15), "mapFileup", "xmlFileup", "imageFileup")


def get_time_stamps(connection):

    cursor = connection.cursor()
    cursor.execute("SELECT * FROM time_stamp")
    time_stamps = [
        Time_stamp(
            record[0],
            record[1],
            record[2],
            record[3],
            record[4]) for record in cursor]
    cursor.close()

    return time_stamps


def update_time_stamps(connection, time_stamps):
    cursor = connection.cursor()

    for time_stamp in time_stamps:
        time_stamp.update_db(cursor)

    connection.commit()
    cursor.close()


def insert_time_stamps(connection, time_stamps):
    cursor = connection.cursor()

    for time_stamp in time_stamps:
        time_stamp.insert_db(cursor)

    connection.commit()
    cursor.close()


def create_time_stamps():
    resultado = []
    resultado.append(
        Time_stamp(
            1,
            (2018,
             2,
             14),
            "mapFile",
            "xmlFile",
            "imageFile"))
    return resultado


def create_time_stamp_table(connection):
    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE time_stamp(\
        emd_entry_id INT PRIMARY KEY,\
        modification DATE NOT NULL,\
        map_file TEXT NOT NULL,\
        xml_file TEXT NOT NULL,\
        image_file TEXT NOT NULL\
        );")
    connection.commit()
    cursor.close()
