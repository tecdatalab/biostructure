'''
Created on 10 mar. 2019

@author: luis98
'''
from classes.Representation import Representation
import pytest
from pytest_postgresql import factories  # Import the ability to use postgresql

representation = Representation(None, None)


def test_0_set_get_id():

    for i in range(4):
        id = i
        representation.set_id(id)
        assert representation.get_id() == id


def test_1_set_get_name():

    for i in range(4):
        name = 'name' + str(i)
        representation.set_name(name)
        assert representation.get_name() == name


def test_3_insert_db(postgresql):

    create_representation_table(postgresql)
    representations = create_representations()
    insert_representations(postgresql, representations)
    result = get_representations(postgresql)
    assert result[0] == Representation(1, "Name 1")


def test_4_update_db(postgresql):

    create_representation_table(postgresql)
    representations = create_representations()
    insert_representations(postgresql, representations)
    representations[0].name = "Name 2"

    update_representations(postgresql, representations)
    result = get_representations(postgresql)
    assert result[0] == Representation(1, "Name 2")


def get_representations(connection):

    cursor = connection.cursor()
    cursor.execute("SELECT * FROM representation")
    representations = [
        Representation(
            record[0],
            record[1]) for record in cursor]
    cursor.close()

    return representations


def update_representations(connection, representations):
    cursor = connection.cursor()

    for representation in representations:
        representation.update_db(cursor)

    connection.commit()
    cursor.close()


def insert_representations(connection, representations):
    cursor = connection.cursor()

    for representation in representations:
        representation.insert_db(cursor)

    connection.commit()
    cursor.close()


def create_representations():
    resultado = []
    resultado.append(Representation(1, "Name 1"))
    return resultado


def create_representation_table(connection):
    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE representation(\
        id SERIAL PRIMARY KEY,\
        name TEXT UNIQUE NOT NULL\
        );")
    connection.commit()
    cursor.close()
