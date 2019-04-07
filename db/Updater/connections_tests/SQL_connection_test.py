'''
Created on 17 mar. 2019

@author: luis98
'''
from classes.Update import Update
from connections.SQL_connection import SQL_connection
import pytest
from pytest_postgresql import factories # Import the ability to use postgresql
import datetime

sql_connection = SQL_connection()

def test_0_commit(postgresql):
    sql_connection.con = postgresql
    cursor_sql = sql_connection.get_cursor()
    create_update_table(cursor_sql)
    sql_connection.commit()
    insert_updates(cursor_sql,create_updates())
    sql_connection.commit()
    cursor_sql.close()
    cursor_sql = sql_connection.get_cursor()
    result = get_updates(cursor_sql)
    assert len(result) == 2

def test_1_is_first_time(postgresql):
    sql_connection.con = postgresql
    cursor_sql = sql_connection.get_cursor()
    create_update_table(cursor_sql)
    sql_connection.commit()
    result = sql_connection.is_first_time()
    assert result == True
    
def test_2_is_first_time(postgresql):
    sql_connection.con = postgresql
    cursor_sql = sql_connection.get_cursor()
    create_update_table(cursor_sql)
    sql_connection.commit()
    insert_updates(cursor_sql,create_updates())
    sql_connection.commit()
    result = sql_connection.is_first_time()
    assert result == False

def test_3_last_update(postgresql):
    sql_connection.con = postgresql
    cursor_sql = sql_connection.get_cursor()
    create_update_table(cursor_sql)
    sql_connection.commit()
    insert_updates(cursor_sql,create_updates())
    sql_connection.commit()
    result = sql_connection.last_update()
    assert result == datetime.datetime(2018,2,14)


def get_updates(cursor):
    cursor.execute("SELECT * FROM update")
    updates = [Update(record[0]) for record in cursor]
    return updates


def insert_updates(cursor, updates):    
    for update in updates:
        update.insert_db(cursor)

def create_updates():
    resultado = []
    resultado.append(Update((2018,2,14)))
    resultado.append(Update((2000,2,14)))
    return resultado

def create_update_table(cursor):
    cursor.execute(
        "CREATE TABLE update(\
        last_update DATE PRIMARY KEY\
        );")

