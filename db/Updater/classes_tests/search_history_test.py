'''
Created on 17 mar. 2019

@author: luis98
'''

from classes.Search_history import Search_history
import pytest
from pytest_postgresql import factories  # Import the ability to use postgresql
from datetime import datetime

search_history = Search_history(
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    None)


def test_0_set_get_id():

    search_history.set_id(1)
    assert search_history.get_id() == 1


def test_1_set_get_date_time():

    search_history.set_date_time(datetime.fromtimestamp(1545730073))
    assert search_history.get_date_time() == datetime.fromtimestamp(1545730073)


def test_2_set_get_ip():

    search_history.set_ip("1.1.1.1")
    assert search_history.get_ip() == "1.1.1.1"


def test_3_set_get_emd_entry_id():

    search_history.set_emd_entry_id(1)
    assert search_history.get_emd_entry_id() == 1


def test_4_set_get_name_file():

    search_history.set_name_file("fileName")
    assert search_history.get_name_file() == "fileName"


def test_5_set_get_contour_level():

    search_history.set_contour_level(1.0)
    assert search_history.get_contour_level() == 1.0


def test_6_set_get_representation_id():

    search_history.set_representation_id(1)
    assert search_history.get_representation_id() == 1


def test_7_set_get_volume_filter_id():

    search_history.set_volume_filter_id(1)
    assert search_history.get_volume_filter_id() == 1


def test_8_set_get_resolution_filter_min():

    search_history.set_resolution_filter_min(1.0)
    assert search_history.get_resolution_filter_min() == 1.0


def test_9_set_get_resolution_filter_max():

    search_history.set_resolution_filter_max(1.0)
    assert search_history.get_resolution_filter_max() == 1.0


def test_10_insert_db(postgresql):
    create_search_history_table(postgresql)
    search_historys = create_search_historys()
    insert_search_historys(postgresql, search_historys)
    result = get_search_historys(postgresql)
    assert result[0] == Search_history(1, datetime.fromtimestamp(
        1545730073), "1.1.1.1", 1, "fileName", 1.0, 1, 1, 1.0, 1.0)


def test_10_update_db(postgresql):
    create_search_history_table(postgresql)
    search_historys = create_search_historys()
    insert_search_historys(postgresql, search_historys)
    search_historys[0].date_time = datetime.fromtimestamp(1545730073)
    search_historys[0].ip = "2.2.3.3"
    search_historys[0].emd_entry_id = 2
    search_historys[0].name_file = "fileNameUp"
    search_historys[0].contour_level = 2.0
    search_historys[0].representation_id = 2
    search_historys[0].volume_filter_id = 2
    search_historys[0].resolution_filter_min = 2.0
    search_historys[0].resolution_filter_max = 2.0

    update_search_historys(postgresql, search_historys)
    result = get_search_historys(postgresql)
    assert result[0] == Search_history(1, datetime.fromtimestamp(
        1545730073), "2.2.3.3", 2, "fileNameUp", 2.0, 2, 2, 2.0, 2.0)


def get_search_historys(connection):

    cursor = connection.cursor()
    cursor.execute("SELECT * FROM search_history")
    search_historys = [
        Search_history(
            record[0],
            record[1],
            record[2],
            record[3],
            record[4],
            record[5],
            record[6],
            record[7],
            record[8],
            record[9]) for record in cursor]
    cursor.close()

    return search_historys


def update_search_historys(connection, search_historys):
    cursor = connection.cursor()

    for search_history in search_historys:
        search_history.update_db(cursor)

    connection.commit()
    cursor.close()


def insert_search_historys(connection, search_historys):
    cursor = connection.cursor()

    for search_history in search_historys:
        search_history.insert_db(cursor)

    connection.commit()
    cursor.close()


def create_search_historys():
    resultado = []
    resultado.append(
        Search_history(1, datetime.fromtimestamp(1545730073),
                       "1.1.1.1", 1, "fileName", 1.0, 1, 1, 1.0, 1.0))
    return resultado


def create_search_history_table(connection):
    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE search_history(\
        id SERIAL PRIMARY KEY,\
        date_time TIMESTAMP NOT NULL,\
        ip TEXT NOT NULL,\
        emd_entry_id INT NOT NULL,\
        name_file TEXT,\
        contour_level FLOAT8,\
        representation_id INT NOT NULL,\
        volume_filter_id INT NOT NULL,\
        resolution_filter_min FLOAT8,\
        resolution_filter_max FLOAT8\
        );")
    connection.commit()
    cursor.close()
