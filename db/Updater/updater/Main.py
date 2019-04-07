from connections.FTP_connection import FTP_connection
from connections.SQL_connection import SQL_connection
from classes.Emd_entry import Emd_entry
from classes.Update import Update
from datetime import date
from classes.Time_stamp import Time_stamp
from psycopg2 import sql
'''
Created on 31 mar. 2019

@author: luis98
'''

def first_time_emd(conec_ftp,conec_sql):
    cursor_sql = conec_sql.get_cursor()
    emds = conec_ftp.get_all_emds_id()
    
    for i in emds:
        temp_emd_entry = Emd_entry()
        temp_emd_entry.create_by_ftp(i,conec_ftp.ftp)
        temp_emd_entry.insert_db(cursor_sql)
        
        temp_time_stamp = Time_stamp (i,date.today(),temp_emd_entry.map_url,temp_emd_entry.xml_url,temp_emd_entry.image_url )
        temp_time_stamp.insert_db(cursor_sql)
    
    update_date = Update(date.today())
    update_date.insert_db(cursor_sql)
    conec_sql.commit()
    cursor_sql.close()
    
def update_emd(conec_ftp,conec_sql):
    cursor_sql = conec_sql.get_cursor()
    
    emds = conec_ftp.get_emds_higher_than_date(conec_sql.last_update())
    
    for i in emds:
        temp_emd_entry = Emd_entry()
        temp_emd_entry.create_by_ftp(i,conec_ftp.ftp)
        temp_time_stamp = Time_stamp (i,date.today(),temp_emd_entry.map_url,temp_emd_entry.xml_url,temp_emd_entry.image_url )
        
        cursor_sql.execute(sql.SQL("SELECT map_information_id FROM emd_entry WHERE id = %s")
        ,[temp_emd_entry.id])
        map_result = [record[0] for record in cursor_sql]
        if len(map_result)==1:
            temp_emd_entry.map_id = map_result[0]
            temp_emd_entry.update_db(cursor_sql)
            temp_time_stamp.update_db(cursor_sql)
        else:
            temp_emd_entry.insert_db(cursor_sql)
            temp_time_stamp.insert_db(cursor_sql)
    update_date = Update(date.today())
    update_date.insert_db(cursor_sql)
    conec_sql.commit()
    cursor_sql.close()
    

def main():
    conec_ftp = FTP_connection()
    conec_ftp.init_connection()
    conec_sql = SQL_connection()
    conec_sql.init_connection()
    
    if conec_sql.is_first_time():
        first_time_emd(conec_ftp, conec_sql)
    else:
        update_emd(conec_ftp,conec_sql)
        
    conec_sql.close_connection()    
    conec_ftp.close_connection()
    print("Finish")
  
if __name__== "__main__":
    main()

