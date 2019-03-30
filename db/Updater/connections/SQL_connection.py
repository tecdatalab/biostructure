'''
Created on 30 mar. 2019

@author: luis98
'''
import psycopg2

class SQL_connection(object):
    '''
    classdocs
    '''
    # Postgres
    __PSQL_HOST = None
    __PSQL_PORT = None
    __PSQL_USER = None
    __PSQL_PASS = None
    __PSQL_DB   = None
    __con = None

    def __init__(self, PSQL_HOST = "localhost", PSQL_PORT = "5432", PSQL_USER ="postgres" , PSQL_PASS ="root", PSQL_DB = "biomolecules_db"):
        self.__PSQL_HOST = PSQL_HOST
        self.__PSQL_PORT = PSQL_PORT
        self.__PSQL_USER = PSQL_USER
        self.__PSQL_PASS = PSQL_PASS
        self.__PSQL_DB = PSQL_DB
        
    def init_connection(self):
        connstr = "host=%s port=%s user=%s password=%s dbname=%s" % (self.__PSQL_HOST, self.__PSQL_PORT, self.__PSQL_USER, self.__PSQL_PASS, self.__PSQL_DB)
        self.__conn = psycopg2.connect(connstr)
        
        
    def commit(self):
        self.__conn.commit()
        
        
    def close_connection(self):
        self.__conn.close()
        
        
    def get_cursor(self):
        return self.__conn.cursor()


    def get_psql_host(self):
        return self.__PSQL_HOST


    def get_psql_port(self):
        return self.__PSQL_PORT


    def get_psql_user(self):
        return self.__PSQL_USER


    def get_psql_pass(self):
        return self.__PSQL_PASS


    def get_psql_db(self):
        return self.__PSQL_DB


    def get_con(self):
        return self.__con


    def set_psql_host(self, value):
        self.__PSQL_HOST = value


    def set_psql_port(self, value):
        self.__PSQL_PORT = value


    def set_psql_user(self, value):
        self.__PSQL_USER = value


    def set_psql_pass(self, value):
        self.__PSQL_PASS = value


    def set_psql_db(self, value):
        self.__PSQL_DB = value


    def set_con(self, value):
        self.__con = value


    def del_psql_host(self):
        del self.__PSQL_HOST


    def del_psql_port(self):
        del self.__PSQL_PORT


    def del_psql_user(self):
        del self.__PSQL_USER


    def del_psql_pass(self):
        del self.__PSQL_PASS


    def del_psql_db(self):
        del self.__PSQL_DB


    def del_con(self):
        del self.__con

    PSQL_HOST = property(get_psql_host, set_psql_host, del_psql_host, "PSQL_HOST's docstring")
    PSQL_PORT = property(get_psql_port, set_psql_port, del_psql_port, "PSQL_PORT's docstring")
    PSQL_USER = property(get_psql_user, set_psql_user, del_psql_user, "PSQL_USER's docstring")
    PSQL_PASS = property(get_psql_pass, set_psql_pass, del_psql_pass, "PSQL_PASS's docstring")
    PSQL_DB = property(get_psql_db, set_psql_db, del_psql_db, "PSQL_DB's docstring")
    con = property(get_con, set_con, del_con, "con's docstring")

'''
test = SQL_connection()
test.init_connection()
cursor = test.get_cursor()
objecto = Volume_filter(1,"name")
objecto.insert_db(cursor)
test.commit()
cursor.close()
test.close_connection()
'''