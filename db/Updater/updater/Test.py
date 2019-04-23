from classes.Descriptor import Descriptor
from connections.SQL_connection import SQL_connection

'''
Created on 8 abr. 2019

@author: luis98
'''

descriptor = Descriptor(1,1,[1,2,3])

conec_sql = SQL_connection()
conec_sql.init_connection()

cursor_sql = conec_sql.get_cursor()
    
descriptor.insert_db(cursor_sql)    
conec_sql.commit()
cursor_sql.close()
conec_sql.close_connection() 
print("Fin")   

