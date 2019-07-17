import sys
sys.path.append('../')
from connections.FTP_connection import FTP_connection
from connections.SQL_connection import SQL_connection
from classes.Emd_entry import Emd_entry
from classes.Update import Update
from datetime import date
from classes.Time_stamp import Time_stamp
from psycopg2 import sql
from generators.Values_generator import get_emd_descriptors, remove_map, download_file
from classes.Descriptor import Descriptor
from generators import Values_generator as valg
from time import time 
from generators import Gif_generator as gifG
import argparse 
'''
Created on 31 mar. 2019

@author: luis98
'''

valg.dir = "../generators/"
gifG.dir = "../generators/"
log_file = None

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

def update_emd(conec_ftp,conec_sql,emd_url_server, emd_url_path, initialEMD, mode, image, descriptor):
    cursor_sql = conec_sql.get_cursor()
    
    if mode == 'c':
        emds = conec_ftp.get_all_emds_id(initialEMD)
    else:
        emds = conec_ftp.get_emds_higher_than_date(conec_sql.last_update(),initialEMD)
    
    print("Total changes",len(emds))
    k=1
    for i in emds:
        temp_emd_entry = Emd_entry()
        temp_emd_entry.emd_url = "http://{0}{1}".format(emd_url_server,emd_url_path)
        temp_emd_entry.create_by_ftp(i,conec_ftp.ftp)
        temp_time_stamp = Time_stamp (i,date.today(),temp_emd_entry.map_url,temp_emd_entry.xml_url,temp_emd_entry.image_url )
        
        cursor_sql.execute(sql.SQL("SELECT map_information_id FROM emd_entry WHERE id = %s")
        ,[temp_emd_entry.id])
        map_result = [record[0] for record in cursor_sql]
        download_file(i)
        
        if image == 'T':
            try:
                imagpath,gifpath = gifG.generateGif(i,temp_emd_entry.map_countour_level)
            except Exception as e:
                log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                log_file.write("Error in the images generation of EMD {0}".format(i))
                log_file.write(str(e))
                log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                print("Error in image generator")
                print(str(e))
        else:
            imagpath,gifpath = None, None
        
        temp_emd_entry.gif_img_3d = gifpath
        temp_emd_entry.png_img_3d = imagpath            
        
        if len(map_result)==1:
            temp_emd_entry.map_id = map_result[0]
            try:
                if image == 'T':
                    temp_emd_entry.update_db(cursor_sql)
                else:
                    temp_emd_entry.update_db_without_images(cursor_sql)
            except Exception as e:
                log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                log_file.write("Error in the update of EMD {0}".format(i))
                log_file.write(str(e))
                log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                print("Error in the update")
                print(str(e))
            
            if descriptor == 'T':
                try:
                    update_descriptor(temp_emd_entry,cursor_sql)
                except Exception as e:
                    log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                    log_file.write("Error in the descriptor generation of EMD {0}".format(i))
                    log_file.write(str(e))
                    log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                    print("Error in descriptor generator")
                    print(str(e))
                    
            temp_time_stamp.update_db(cursor_sql)
        else:
            try:
                temp_emd_entry.insert_db(cursor_sql)
            except Exception as e:
                log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                log_file.write("Error in the insert of EMD {0}".format(i))
                log_file.write(str(e))
                log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                print("Error in the insert")
                print(str(e))
            
            if descriptor == 'T':
                try:
                    insert_descriptor(temp_emd_entry,cursor_sql)
                except Exception as e:
                    log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                    log_file.write("Error in the descriptor generation of EMD {0}".format(i))
                    log_file.write(str(e))
                    log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
                    print("Error in descriptor generator")
                    print(str(e))
                    
            temp_time_stamp.insert_db(cursor_sql)
        
        remove_map(i)
        conec_sql.commit()
        print("Actual execution : {0} with EMD: {1}".format(k,i))
        k+=1
    
    if(conec_sql.last_update().date()!=date.today()):
        update_date = Update(date.today())
        update_date.insert_db(cursor_sql)
    
    conec_sql.commit()
    cursor_sql.close()
    

def main(emd_url_server, emd_url_path, initialEMD, mode, image, descriptor):
    valg.emd_url = "http://{0}{1}".format(emd_url_server,emd_url_path)
    conec_ftp = FTP_connection(ftp_server=emd_url_server)
    conec_ftp.init_connection(initial_directory="{0}{1}".format(emd_url_path,"/structures/"))
    conec_sql = SQL_connection()
    conec_sql.init_connection()
    
    update_emd(conec_ftp,conec_sql,emd_url_server, emd_url_path, initialEMD, mode, image, descriptor)
        
    conec_sql.close_connection()    
    conec_ftp.close_connection()
    print("Finish")
  
if __name__== "__main__":
    parser = argparse.ArgumentParser(description='Calculation of data for EMD.')
    parser.add_argument('-l','--log', help='Log file name.', default='log.txt')
    parser.add_argument('-ie','--initialEMD', help='Initial EMD for execution.', default='0001')
    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument('-m', '--mode', choices=['c','u'], help='Execution mode, where "c" is complete and "u" update.', required=True)
    required_arguments.add_argument('-i', '--image', choices=['Y','N'], help='Generation of images and gif where "Y" is yes and "N" is no.', required=True)
    required_arguments.add_argument('-d', '--descriptor', choices=['Y','N'], help='Generation of descriptor where "Y" is yes and "N" is no.', required=True)
    args = parser.parse_args()
    
    log_file = open(args.log, 'a+')  # open file in append mode
    log_file.write("========================================")
    ini = time() 
    main("ftp.wwpdb.org", "/pub/emdb", args.initialEMD, args.mode, args.image, args.descriptor)
    final = time() 
    ejec = final - ini
    log_file.close()
    print ('Execution time:', ejec)
