import sys, argparse, os
sys.path.append('../')
import datetime

from connections.ftp_connection import FTP_connection
from connections.sql_connection import SQL_connection
from classes.update import Update
from classes.binnacle import Binnacle
from main_pdb import main

def main_script_pdb(connec_sql, connec_ftp, cathComplex, cathChain, cathDomain, atomic, startPdb, amount, atomicEmd):
    cursor_sql = connec_sql.get_cursor()

    update_temp = Update(datetime.datetime.now())
    update_temp.insert_update_db(connec_sql, cursor_sql)
    connec_sql.commit()

    temp_binnacle = Binnacle()
    temp_binnacle.set_last_update(cursor_sql)
    temp_pdb_id = temp_binnacle.get_last_pdb_id_update(cursor_sql)
    
    connec_sql.close_connection()
    connec_ftp.close_connection()

    if(temp_pdb_id == None):
        main(cathComplex, cathChain, cathDomain, atomic, startPdb,  startPdb + amount, atomicEmd)
    else:
        index_pdb = connec_ftp.get_pdb_by_index(temp_pdb_id)
        main(cathComplex, cathChain, cathDomain, atomic, index_pdb,  index_pdb + amount, atomicEmd)

def init_script_pdb(cathComplex, cathChain, cathDomain, atomic, startPdb, amount, atomicEmd):
    connec_sql = SQL_connection()
    connec_sql.init_connection()
    
    connec_ftp = FTP_connection()
    connec_ftp.init_connection()
    
    main_script_pdb(connec_sql, connec_ftp, cathComplex, cathChain, cathDomain, atomic, int(startPdb), int(amount), atomicEmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculation of data for PDB.')
    parser.add_argument(
        '-l',
        '--log',
        help='Log file name.',
        default='log.txt')

    parser.add_argument(
        '-ho',
        '--hour',
        help='Time schedule.',
        default='1')

    parser.add_argument(
        '-ia',
        '--initialAtomic',
        help='Initial Atomic for execution.',
        default='0')

    parser.add_argument(
        '-am',
        '--amount',
        help='Amount of pdb to execute.',

        default='1')

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

    
    args = parser.parse_args()
    init_script_pdb(args.cathComplex, args.cathChain, args.cathDomain, args.atomic, args.initialAtomic, args.amount, args.atomicEmd)