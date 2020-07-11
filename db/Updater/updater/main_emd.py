import sys, argparse, os
sys.path.append('../')
from time import time
from generators import gif_generator as gifG
from generators import values_generator as valg
from classes.descriptor import Descriptor
from generators.values_generator import get_emd_descriptors, remove_map, download_file
from classes.time_stamp import Time_stamp
from datetime import date
from classes.update import Update
from classes.binnacle import Binnacle
from classes.emd_entry import Emd_entry
from connections.sql_connection import SQL_connection
from connections.ftp_connection import FTP_connection


from constants.constants import *
from utilities.utility import *
from utilities.log import Log
'''
Created on 31 mar. 2019
@author: luis98

Last modified by @dnnxl
Date on 24 may 2020
'''

valg.dir = "../generators/"
gifG.dir = "../generators/"
log_file = None
csv_file = 'emd.csv'
attempt_emd = 5

def insert_update_descriptors(emd_entry_p, cursor_sql):
    result = get_emd_descriptors(
        emd_entry_p.id,
        emd_entry_p.map_countour_level,
        emd_entry_p.map_std)
    
    for i in range(len(result)):
        des_temp = Descriptor(emd_entry_p.id, i + 1, result[i])
        des_temp.insert_update_db(cursor_sql)

def update_emd(connec_ftp,connec_sql,initialEMD,mode,image,descriptor,finalEMD):
    
    cursor_sql = connec_sql.get_cursor()
    if mode == 'c':
        emds = connec_ftp.get_all_emds_id(initialEMD, finalEMD)
    else:
        emds = connec_ftp.get_emds_higher_than_date(connec_sql.last_update(), initialEMD, finalEMD)

    temp_binnacle = Binnacle()
    temp_binnacle.set_last_update(cursor_sql)
    
    print("Total changes", len(emds))

    k = 1
    for i in emds:
        
        temp_binnacle.set_emd_id(int(i))
        temp_emd_attempt = temp_binnacle.get_attempt_emd(cursor_sql)
        
        if(temp_emd_attempt == None):
            temp_binnacle.insert_binnacle_emd(cursor_sql)
            connec_sql.commit()
        elif(temp_emd_attempt > attempt_emd):
            continue
        else:
            temp_binnacle.update_binnacle_emd(cursor_sql)

        ini_time = time()
        ini_memory = memory()

        print("Actual execution : {0} with EMD: {1}".format(k, i))
        temp_emd_entry = Emd_entry()
        temp_emd_entry.create_by_ftp(i, connec_ftp.get_ftp())
        temp_time_stamp = Time_stamp(i, date.today(), temp_emd_entry.map_url, temp_emd_entry.xml_url, temp_emd_entry.image_url)
        
        if image == 'Y' or descriptor == 'Y':
            print("Download file")
            download_file(i)

        if image == 'Y':
            try:
                temp_emd_entry.create_gif()

            except Exception as e:
                log_file.generate_error_message(TypeErrorEmd.EDG.value.format(i), TypeErrorEmd.EDG.name, e)

            try:
                temp_emd_entry.insert_update_db(cursor_sql)

            except Exception as e:
                log_file.generate_error_message(TypeErrorEmd.EIN.value.format(i), TypeErrorEmd.EIN.name, e)
        else:
            try:
                temp_emd_entry.insert_without_images_db(cursor_sql)

            except Exception as e:
                log_file.generate_error_message(TypeErrorEmd.EIN.value.format(i), TypeErrorEmd.EIN.name, e)

        if descriptor == 'Y':
            if temp_emd_entry.map_countour_level is not None and temp_emd_entry.map_std is not None:
                try:
                    insert_update_descriptors(temp_emd_entry, cursor_sql)

                except Exception as e:
                    log_file.generate_error_message(TypeErrorEmd.EID.value.format(i), TypeErrorEmd.EID.name, e)

            else:
                log_file.generate_error_message(TypeErrorEmd.ECS.value.format(i), TypeErrorEmd.ECS.name, None)

        try:
            temp_time_stamp.insert_update_db(cursor_sql)
        except Exception as e:
            log_file.generate_error_message(TypeErrorEmd.ETS.value.format(i), TypeErrorEmd.ETS.name, e)

        end_time = time()
        end_memory = memory()
        ejec_memory = end_memory - ini_memory
        ejec_time = end_time - ini_time
        row_content = [i, ejec_memory, ejec_time]
        append_list_as_row(csv_file, row_content)
        
        remove_map(i)
        connec_sql.commit()
        k += 1
        
    update_temp = Update(date.today())
    update_temp.insert_update_db(connec_sql, cursor_sql)
    connec_sql.commit()
    cursor_sql.close()

def main(initialEMD, mode, image, descriptor, finalEMD):
    connec_ftp = FTP_connection()
    connec_ftp.init_connection()
    connec_sql = SQL_connection()
    connec_sql.init_connection()

    update_emd(connec_ftp, connec_sql, initialEMD, mode, image, descriptor, finalEMD)

    connec_sql.close_connection()
    connec_ftp.close_connection()
    print("Finish")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculation of data for EMD.')
    parser.add_argument(
        '-l',
        '--log',
        help='Log file name.',
        default='updater.log')
    parser.add_argument(
        '-ie',
        '--initialEMD',
        help='Initial EMD for execution.',
        default='0001')
    parser.add_argument(
        '-fe',
        '--finalEMD',
        help='Final EMD for execution (use "inf" for a complete execution).',
        default='inf')
    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument(
        '-m',
        '--mode',
        choices=[
            'c',
            'u'],
        help='Execution mode, where "c" is complete and "u" update.',
        required=True)
    required_arguments.add_argument(
        '-i',
        '--image',
        choices=[
            'Y',
            'N'],
        help='Generation of images and gif where "Y" is yes and "N" is no.',
        required=True)
    required_arguments.add_argument(
        '-d',
        '--descriptor',
        choices=[
            'Y',
            'N'],
        help='Generation of descriptor where "Y" is yes and "N" is no.',
        required=True)
    args = parser.parse_args()

    log_file = Log(args.log)  
    log_file.init_log_file(__file__)
    log_file.generate_info_message(TypeMessage.MS1.value, TypeMessage.MS1.name)
    ini = time() # The initial execution time 
    main(args.initialEMD, args.mode, args.image, args.descriptor, args.finalEMD)
    final = time()
    ejec = final - ini
    print('Execution time:', ejec)
    log_file.generate_info_message(TypeMessage.MS2.value, TypeMessage.MS2.name)
