import sys, argparse, os
import datetime
sys.path.append('../')
from time import time
from connections.sql_connection import SQL_connection
from connections.ftp_connection import FTP_connection
from classes.atomic_structure import Atomic_structure
from classes.binnacle import Binnacle
from datetime import date
from classes.update import Update

from utilities.log import Log
from constants.constants import *
'''
Created on 31 mar. 2019
@author: luis98

Last modified: 16 may 2020
@By: dnnxl
'''

log_file = None
attempt_pdb = 5

def update_pdb(connec_ftp, connec_sql, cathComplex, cathChain, cathDomain, atomic, atomicEmd, initialAtomic, finalAtomic):
    
    cursor_sql = connec_sql.get_cursor()

    update_temp = Update(datetime.datetime.now())
    update_temp.insert_update_db(connec_sql, cursor_sql)
    connec_sql.commit()

    temp_binnacle = Binnacle()
    temp_binnacle.set_last_update(cursor_sql)

    #Create atomic structures
    if atomic == 'Y':
        # Complete mode
        if(int(initialAtomic) == 0 and float(finalAtomic) == float("inf")):
            all_pdb = connec_ftp.get_all_pdb()
        else:
            all_pdb = connec_ftp.get_range_pdb(int(initialAtomic), float(finalAtomic))

        print("Total changes atomic structures", len(all_pdb))
        domain_dic = connec_ftp.get_all_cath_domain_boundarie_dic()
        k = 1
        for i in all_pdb:

            temp_binnacle.set_pdb_id(i.id_code)
            temp_pdb_attempt = temp_binnacle.get_attempt_pdb(cursor_sql)
            if(temp_pdb_attempt == None):
                temp_binnacle.insert_binnacle_pdb(cursor_sql)
                connec_sql.commit()
            elif(temp_pdb_attempt > attempt_pdb):
                continue
            else:
                temp_binnacle.update_binnacle_pdb(cursor_sql)
                connec_sql.commit()

            print("---------------------------------------------------------")
            print ("Actual execution {0} with PDB {1}".format(k, i.id_code))
            try:
                temp = Atomic_structure()
                temp.create_by_online_file(i.id_code, domain_dic.get(i.id_code))                    
                temp.insert_update_db_complex(cursor_sql)
            except Exception as e:
                log_file.generate_error_message(TypeErrorPdb.EEP.value.format(k, i.id_code), TypeErrorPdb.EEP.name, e)
            k += 1
            connec_sql.commit()

    #Create cathComplex
    if cathComplex == 'Y':
        all_cathChain = connec_ftp.get_all_cath_complex()
        print("Total changes cath complex atomic structures", len(all_cathChain))
        k = 1
        for i in all_cathChain:
            print("---------------------------------------------------------")
            print ("Actual cath complex {0} with complex {1}".format(k, i.id_code))
            try:
                i.insert_update_db(cursor_sql)
            except Exception as e:
                log_file.generate_error_message(TypeErrorPdb.ECC.value.format(k, i.id_code), TypeErrorPdb.ECC.name, e)
            k += 1
            connec_sql.commit()
    
    #Create cathChain
    if cathChain == 'Y':
        all_cathChain = connec_ftp.get_all_cath_chain()
        print("Total changes cath chain atomic structures", len(all_cathChain))
        k = 1
        for i in all_cathChain:
            print("---------------------------------------------------------")
            print ("Actual cath chain {0} with chain {1}".format(k, i.id_code))
            try:
                i.insert_update_db(cursor_sql)
            except Exception as e:
                log_file.generate_error_message(TypeErrorPdb.ECH.value.format(k, i.id_code), TypeErrorPdb.ECH.name, e)
            k += 1
            connec_sql.commit()
    
    #Create cathDomain
    if cathDomain == 'Y':
        all_cathDomain = connec_ftp.get_all_cath_domain()
        print("Total changes cath domain atomic structures", len(all_cathDomain))
        k = 1
        for i in all_cathDomain:
            print("---------------------------------------------------------")
            print ("Actual cath domain {0} with domain {1}".format(k, i.id_code))
            try:
                i.insert_update_db(cursor_sql)
            except Exception as e:
                pass
                log_file.generate_error_message(TypeErrorPdb.ECD.value.format(k, i.id_code), TypeErrorPdb.ECD.name, e)
            k += 1
            connec_sql.commit()
    
    #Create atomic x emd
    if atomicEmd == 'Y':
        all_se = connec_ftp.get_all_structure_x_emd_entry()
        print("Total changes atomic structure x emd entry", len(all_se))
        k = 1
        for i in all_se:
            print("---------------------------------------------------------")
            print ("Actual atomic structure x emd_entry {0} with atomic_structure {1}, emd_entry {2}".format(k, i.atomic_structure , i.emd_entry))
            try:
                i.insert_update_db(cursor_sql)
            except Exception as e:
                log_file.generate_error_message(TypeErrorPdb.EAS.value.format(k, i.atomic_structure , i.emd_entry), TypeErrorPdb.EAS.name, e)
            k += 1
            connec_sql.commit()
    
    connec_sql.commit()
    cursor_sql.close()

def main(cathComplex, cathChain, cathDomain, atomic, initialAtomic, finalAtomic, atomicEmd):
    
    connec_ftp = FTP_connection()
    connec_ftp.init_connection()
    connec_sql = SQL_connection()
    connec_sql.init_connection()

    update_pdb(connec_ftp, connec_sql, cathComplex, cathChain, cathDomain, atomic, atomicEmd, initialAtomic, finalAtomic)

    connec_sql.close_connection()
    connec_ftp.close_connection()
    print("Finish")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculation of data for PDB.')
    parser.add_argument(
        '-l',
        '--log',
        help='Log file name.',
        default='log.txt')

    parser.add_argument(
        '-ia',
        '--initialAtomic',
        help='Initial Atomic for execution.',
        default='0')
    parser.add_argument(
        '-fa',
        '--finalAtomic',
        help='Final Atomic for execution (use "inf" for a complete execution).',
        default='inf')

    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument(
        '-cco',
        '--cathComplex',
        choices=[
            'Y',
            'N'],
        help='Generate cath atomic structures for complex.',
        required=True)
    required_arguments.add_argument(
        '-cc',
        '--cathChain',
        choices=[
            'Y',
            'N'],
        help='Generate cath atomic structures for chain.',
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
            'N'],
        help='Generation atomic structures.',
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

    log_file = Log(args.log)  # open file in append mode
    log_file.init_log_file(__file__)
    log_file.generate_info_message(TypeMessage.MS1.value, TypeMessage.MS1.name)
    ini = time()
    main(args.cathComplex, args.cathChain, args.cathDomain, args.atomic, args.initialAtomic, args.finalAtomic, args.atomicEmd)
    final = time()
    ejec = final - ini
    print('Execution time:', ejec)
    log_file.generate_info_message(TypeMessage.MS2.value, TypeMessage.MS2.name)
