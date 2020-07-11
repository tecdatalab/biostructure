'''
@author dnnxl
Created 28 may 2020
'''
import sys, argparse, os
sys.path.append('../')
import datetime

from crontab import CronTab
from classes.update import Update
from connections.sql_connection import SQL_connection

directory = os.path.abspath(os.getcwd())

def initial():
    cron = CronTab(user=True)
    cron.remove_all()
    cron.write()

    connec_sql = SQL_connection()
    connec_sql.init_connection()
    cursor_sql = connec_sql.get_cursor()

    update_temp = Update(datetime.datetime.now())
    update_temp.insert_last_update_db(connec_sql, cursor_sql)
    connec_sql.commit()
    connec_sql.close_connection()

def schedule(cathComplex, cathChain, cathDomain, atomic, initialAtomic, amount, atomicEmd, hour):
    cron = CronTab(user=True)
    temp_command = "sudo python3 {}/script_pdb.py -cco {} -cc {} -cd {} -a {} -ae {} -am {} -ia {}".format(directory, cathComplex, cathChain, cathDomain, atomic, atomicEmd, amount, initialAtomic)
    job = cron.new(command=temp_command)
    job.every(hour).hours()
    cron.write()

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
        '-am',
        '--amount',
        help='Amount of pdb to execute.',

        default='1')
    
    parser.add_argument(
        '-ho',
        '--hour',
        help='Time schedule.',
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
    print(args)
    initial()

    initial()
    schedule(
        args.cathComplex, 
        args.cathChain, 
        args.cathDomain, 
        args.atomic, 
        args.initialAtomic, 
        args.amount, 
        args.atomicEmd, 
        args.hour)