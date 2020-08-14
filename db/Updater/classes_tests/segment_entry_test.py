'''
Created on 7 aug. 2020

@author: dnnxl
'''

from classes.segment_entry import Segment_entry, List_segments
import pytest
from pytest_postgresql import factories

segment_entry = Segment_entry(None, None, 
                              None, None, 
                              None, None)

list_segments_entry = List_segments(None)

def test_0_set_get_id():

    id = 1
    segment_entry.segment_id = id
    assert segment_entry.segment_id == id

def test_1_set_get_map_id():

    id_map = 1
    segment_entry.segment_map_id = id_map
    assert segment_entry.segment_map_id == id_map

def test_2_set_get_algorithm_id():

    id_algorithm = 1
    segment_entry.segment_id_algorithm = id_algorithm
    assert segment_entry.segment_id_algorithm == id_algorithm

def test_3_set_get_volume():

    volume = 1000
    segment_entry.segment_volume = volume
    assert segment_entry.segment_volume == volume

def test_4_set_get_contour_level():

    contour_level = 1
    segment_entry.segment_contour_level = contour_level
    assert segment_entry.segment_contour_level == contour_level

def test_5_set_get_path():

    path = 'path'
    segment_entry.segment_path = path
    assert segment_entry.segment_path == path

def insert_segments_entries(connection, segment_entries):
    cursor = connection.cursor()

    for segment_entry in segment_entries:
        segment_entry.insert_db(cursor)

    connection.commit()
    cursor.close()

def create_segment_entries():
    result = []
    result.append(
        Segment_entry(1, 1, 1, 123, 52, 'path'))
    return result

def get_segment_entries(connection):

    cursor = connection.cursor()
    cursor.execute(
        "SELECT * FROM segment_entry")
    segment_entries = [
        Segment_entry(
            record[0],
            record[1],
            record[2],
            record[3],
            record[4],
            record[5]) for record in cursor]
    cursor.close()

    return segment_entries

def test_6_insert_db(connection):
    segment_entries = create_segment_entries()
    insert_segments_entries(connection, segment_entries)
    result = get_segment_entries(connection)
    assert result[0] == Segment_entry(1, 1, 1, 123, 52, 'path')

def test_7_set_get_list_segment():
    segments = [Segment_entry(1, 1, 1, 123, 52, 'path')]
    list_segments_entry.list_segments = segments
    assert list_segments_entry.list_segments == segments