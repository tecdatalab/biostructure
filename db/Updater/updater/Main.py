from connections.FTP_connection import FTP_connection
from connections.SQL_connection import SQL_connection
from classes.Emd_entry import Emd_entry
'''
Created on 31 mar. 2019

@author: luis98
'''

def main():
    conec_ftp = FTP_connection()
    conec_ftp.init_connection()
    conec_sql = SQL_connection()
    conec_sql.init_connection()
    cursor_sql = conec_sql.get_cursor()
    emds = conec_ftp.get_all_emds_id()
    
    for i in emds:
        temp = Emd_entry()
        temp.create_by_ftp(i,conec_ftp.ftp)
        temp.insert_db(cursor_sql)
    
    conec_sql.commit()
    cursor_sql.close()
    conec_sql.close_connection()    
    conec_ftp.close_connection()
  
if __name__== "__main__":
    main()

