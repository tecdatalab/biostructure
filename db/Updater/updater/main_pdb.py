import sys, argparse, os
sys.path.append('../')
from time import time
from classes.update import Update
from connections.sql_connection import SQL_connection
from connections.ftp_connection import FTP_connection
from classes.atomic_structure import Atomic_structure

'''
Created on 31 mar. 2019
@author: luis98
'''

log_file = None

def update_log_file():
    log_file.flush()
    os.fsync(log_file.fileno())

def generate_error_message(personalized_error, compute_error):
    log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n")
    log_file.write(personalized_error + "\n")
    log_file.write(str(compute_error) + "\n")
    log_file.write("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n")
    print(personalized_error+ "\n")
    print(str(compute_error) + "\n")
    update_log_file()

def update_pdb(connec_ftp, connec_sql, cathComplex, cathDomain, atomic, atomicEmd):
    
    cursor_sql = connec_sql.get_cursor()
    if atomic == 'Y':
        all_pdb = connec_ftp.get_all_pdb()
        print("Total changes", len(all_pdb))
        domain_dic = connec_ftp.get_all_cath_domain_boundarie_dic()
        k = 1
        for i in all_pdb:
            try:
                temp = Atomic_structure()
                temp.create_by_online_file(i.data1, domain_dic.get(i.data1))
                temp.insert_update_db_without_descriptors(cursor_sql)
            except Exception as e:
                generate_error_message("Error in the generation of PDB {0}".format(i),e)
            print ("Actual execution {0} with PDB {1}".format(i, k))
            k += 1

            connec_sql.commit()

    elif atomic == 'YD':
        all_pdb = connec_ftp.get_all_pdb()
        print("Total changes", len(all_pdb))
        domain_dic = connec_ftp.get_all_cath_domain_boundarie_dic()
        k = 1
        for i in all_pdb:
            try:
                temp = Atomic_structure()
                temp.create_by_online_file(i.data1, domain_dic.get(i.data1))
                temp.calculate_descriptors()
                temp.insert_update_db(cursor_sql)
            except Exception as e:
                generate_error_message("Error in the generation of PDB {0}".format(i),e)
            print ("Actual execution {0} with PDB {1}".format(i, k))
            k += 1

            connec_sql.commit()

    connec_sql.commit()
    cursor_sql.close()

def main(cathComplex, cathDomain, atomic, atomicEmd):
    
    connec_ftp = FTP_connection()
    connec_ftp.init_connection()
    connec_sql = SQL_connection()
    connec_sql.init_connection()

    update_pdb(connec_ftp, connec_sql, cathComplex, cathDomain, atomic, atomicEmd)

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
    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument(
        '-cc',
        '--cathComplex',
        choices=[
            'Y',
            'N'],
        help='Generate cath atomic structures for complex.',
        required=True)
    required_arguments.add_argument(
        '-cd',
        '--cathDomain',
        choices=[
            'Y',
            'N'],
        help='Generate cath atomic structures for domain.',
        required=True)
    required_arguments.add_argument(
        '-a',
        '--atomic',
        choices=[
            'Y',
            'N',
            'YD',
            'ND'],
        help='Generation atomic structures. Letter D means add descriptors',
        required=True)
    required_arguments.add_argument(
        '-ae',
        '--atomicEmd',
        choices=[
            'Y',
            'N'],
        help='Generation connection with pdb and emd.',
        required=True)
    args = parser.parse_args()

    log_file = open(args.log, 'a+')  # open file in append mode
    log_file.write("====================Start====================\n")
    update_log_file()
    ini = time()
    main(args.cathComplex, args.cathDomain, args.atomic, args.atomicEmd)
    final = time()
    ejec = final - ini
    log_file.close()
    print('Execution time:', ejec)