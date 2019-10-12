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
from classes.emd_entry import Emd_entry
from connections.sql_connection import SQL_connection
from connections.ftp_connection import FTP_connection

'''
Created on 31 mar. 2019
@author: luis98
'''

valg.dir = "../generators/"
gifG.dir = "../generators/"
log_file = None

def update_log_file():
    log_file.flush()
    os.fsync(log_file.fileno())

def insert_update_descriptors(emd_entry_p, cursor_sql):
    result = get_emd_descriptors(
        emd_entry_p.id,
        emd_entry_p.map_countour_level,
        emd_entry_p.map_std)
    
    for i in range(len(result)):
        des_temp = Descriptor(emd_entry_p.id, i + 1, result[i])
        des_temp.insert_update_db(cursor_sql)

def generate_error_message(personalized_error, compute_error):
    log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n")
    log_file.write(personalized_error + "\n")
    log_file.write(str(compute_error) + "\n")
    log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n")
    print(personalized_error+ "\n")
    print(str(compute_error) + "\n")
    update_log_file()

def update_emd(connec_ftp,connec_sql,initialEMD,mode,image,descriptor,finalEMD):
    
    cursor_sql = connec_sql.get_cursor()
    if mode == 'c':
        emds = connec_ftp.get_all_emds_id(initialEMD, finalEMD)
    else:
        emds = connec_ftp.get_emds_higher_than_date(connec_sql.last_update(), initialEMD, finalEMD)

    print("Total changes", len(emds))

    k = 1
    for i in emds:
        print("Actual execution : {0} with EMD: {1}".format(k, i))
        temp_emd_entry = Emd_entry()
        temp_emd_entry.create_by_ftp(i, connec_ftp.ftp)
        temp_time_stamp = Time_stamp(i, date.today(), temp_emd_entry.map_url, temp_emd_entry.xml_url, temp_emd_entry.image_url)
        
        if image == 'Y' or descriptor == 'Y':
            download_file(i)

        if image == 'Y':
            try:
                temp_emd_entry.create_gif()
            except Exception as e:
                generate_error_message("Error in the images generation of EMD {0}".format(i),e)
                temp_emd_entry.insert_update_db(cursor_sql)
        else:
            temp_emd_entry.insert_without_images_db(cursor_sql)
        
        if descriptor == 'Y':
            if temp_emd_entry.map_countour_level is not None and temp_emd_entry.map_std is not None:
                try:
                    insert_update_descriptors(temp_emd_entry, cursor_sql)
                except Exception as e:
                    generate_error_message("Error in the descriptor generation of EMD {0}".format(i),e)
            else:
                generate_error_message("Error in the descriptor generation of EMD {0}.Countour level or std not exist.".format(i),None)
        
        temp_time_stamp.insert_update_db(cursor_sql)
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
        default='log.txt')
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

    log_file = open(args.log, 'a+')  # open file in append mode
    log_file.write("====================Start====================\n")
    update_log_file()
    ini = time()
    main(args.initialEMD, args.mode, args.image, args.descriptor, args.finalEMD)
    final = time()
    ejec = final - ini
    log_file.close()
    print('Execution time:', ejec)