from connections.FTP_connection import FTP_connection
from connections.SQL_connection import SQL_connection
from classes.Emd_entry import Emd_entry
from classes.Update import Update
from datetime import date
from classes.Time_stamp import Time_stamp
from psycopg2 import sql
from generators.Values_generator import get_emd_descriptors, remove_map, download_file
from classes.Type_descriptor import Type_descriptor
from classes.Descriptor import Descriptor
from generators import Values_generator as valg
from time import time 
from generators import Gif_generator as gifG
 
'''
Created on 31 mar. 2019

@author: luis98
'''

valg.dir = "../generators/"
gifG.dir = "../generators/"

def create_type_descriptors(cursor_sql):
    
    type = Type_descriptor(1, "EMDB Contour", "EMDB Contour")
    type.insert_db(cursor_sql)
    
    type = Type_descriptor(2, "EMDB Contour + 1/3 core", "EMDB Contour + 1/3 core")
    type.insert_db(cursor_sql)
    
    type = Type_descriptor(3, "EMDB Contour + 2/3 core", "EMDB Contour + 2/3 core")
    type.insert_db(cursor_sql)
    
    type = Type_descriptor(4, "EMDB Contour + 1/3 + 2/3 core", "EMDB Contour + 1/3 + 2/3 core")
    type.insert_db(cursor_sql)
    
    type = Type_descriptor(5, "EMDB Contour + 1 std dev", "EMDB Contour + 1 std dev")
    type.insert_db(cursor_sql)
    
def insert_descriptor(emd_entry_p,cursor_sql):
    result = get_emd_descriptors(emd_entry_p.id,emd_entry_p.map_countour_level,emd_entry_p.map_std)
    
    for i in range(len(result)):
        des_temp = Descriptor(emd_entry_p.id,i+1,result[i])
        des_temp.insert_db(cursor_sql)

def update_descriptor(emd_entry_p,cursor_sql):
    result = get_emd_descriptors(emd_entry_p.id,emd_entry_p.map_countour_level,emd_entry_p.map_std)
    
    for i in range(len(result)):
        des_temp = Descriptor(emd_entry_p.id,i+1,result[i])
        des_temp.update_db(cursor_sql)


def first_time_emd(conec_ftp,conec_sql,emd_url_server, emd_url_path):
    emds = conec_ftp.get_all_emds_id()
    cursor_sql = conec_sql.get_cursor()
    
    create_type_descriptors(cursor_sql)
    print("Total inserts",len(emds))
    for i in emds:
        temp_emd_entry = Emd_entry()
        temp_emd_entry.emd_url = "http://{0}{1}".format(emd_url_server,emd_url_path)
        temp_emd_entry.create_by_ftp(i,conec_ftp.ftp)
        temp_emd_entry.insert_db(cursor_sql)
        
        download_file(i)
        gifG.generateGif(i)
        insert_descriptor(temp_emd_entry,cursor_sql)
        remove_map(i)
        
        temp_time_stamp = Time_stamp (i,date.today(),temp_emd_entry.map_url,temp_emd_entry.xml_url,temp_emd_entry.image_url )
        temp_time_stamp.insert_db(cursor_sql)
        conec_sql.commit()
        print(i)
    
    update_date = Update(date.today())
    update_date.insert_db(cursor_sql)
    conec_sql.commit()
    cursor_sql.close()
    
def update_emd(conec_ftp,conec_sql,emd_url_server, emd_url_path):
    cursor_sql = conec_sql.get_cursor()
    
    emds = conec_ftp.get_emds_higher_than_date(conec_sql.last_update())
    print("Total changes",len(emds))
    for i in emds:
        temp_emd_entry = Emd_entry()
        temp_emd_entry.emd_url = "http://{0}{1}".format(emd_url_server,emd_url_path)
        temp_emd_entry.create_by_ftp(i,conec_ftp.ftp)
        temp_time_stamp = Time_stamp (i,date.today(),temp_emd_entry.map_url,temp_emd_entry.xml_url,temp_emd_entry.image_url )
        
        cursor_sql.execute(sql.SQL("SELECT map_information_id FROM emd_entry WHERE id = %s")
        ,[temp_emd_entry.id])
        map_result = [record[0] for record in cursor_sql]
        if len(map_result)==1:
            temp_emd_entry.map_id = map_result[0]
            temp_emd_entry.update_db(cursor_sql)
            temp_time_stamp.update_db(cursor_sql)
            download_file(i)
            gifG.generateGif(i)
            update_descriptor(temp_emd_entry,cursor_sql)
            remove_map(i)
        else:
            temp_emd_entry.insert_db(cursor_sql)
            temp_time_stamp.insert_db(cursor_sql)
            download_file(i)
            gifG.generateGif(i)
            insert_descriptor(temp_emd_entry,cursor_sql)
            remove_map(i)
        
        conec_sql.commit()
        print(i)
    
    if(conec_sql.last_update().date()!=date.today()):
        update_date = Update(date.today())
        update_date.insert_db(cursor_sql)
    
    conec_sql.commit()
    cursor_sql.close()
    

def main(emd_url_server, emd_url_path):
    valg.emd_url = "http://{0}{1}".format(emd_url_server,emd_url_path)
    conec_ftp = FTP_connection(ftp_server=emd_url_server)
    conec_ftp.init_connection(initial_directory="{0}{1}".format(emd_url_path,"/structures/"))
    conec_sql = SQL_connection()
    conec_sql.init_connection()
    
    if conec_sql.is_first_time():
        first_time_emd(conec_ftp,conec_sql,emd_url_server, emd_url_path)
    else:
        update_emd(conec_ftp,conec_sql,emd_url_server, emd_url_path)
        
    conec_sql.close_connection()    
    conec_ftp.close_connection()
    print("Finish")
  
if __name__== "__main__":
    ini = time() 
    main("ftp.pdbj.org","/pub/emdb")
    final = time() 
    ejec = final - ini
    print ('Execution time:', ejec)


