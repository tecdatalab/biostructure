'''
@author dnnxl
Created 28 may 2020

Create a task using crontab for program
pdb. It runs for each hour, depending the 
parameter hour.

'''
import sys, argparse, os
sys.path.append('../')
import datetime

from crontab import CronTab
from classes.update import Update
from connections.sql_connection import SQL_connection

directory = os.path.abspath(os.getcwd())
crontab_id = 'script_pdb'

def initial():
    cron = CronTab(user=True)
    cron.remove_all(comment=crontab_id)
    cron.write()

    connec_sql = SQL_connection()
    connec_sql.init_connection()
    cursor_sql = connec_sql.get_cursor()
    
    # Format YYYY-MM-DD HH:MM:SS
    update_temp = Update(datetime.datetime.now())
    update_temp.insert_last_update_db(connec_sql, cursor_sql)
    connec_sql.commit()
    connec_sql.close_connection()

def schedule(log, attemptPdb, cathComplex, cathChain, cathDomain, atomic, initialAtomic, amount, atomicEmd, hour):
    cron = CronTab(user=True)
    temp_command = "sudo python3 {}/script_pdb.py -l {} -cco {} -cc {} -cd {} -a {} -ae {} -am {} -ia {} -at {}".format(directory, log, cathComplex, cathChain, cathDomain, atomic, atomicEmd, amount, initialAtomic, attemptPdb)
    job = cron.new(command=temp_command, comment=crontab_id)
    job.every(hour).hours()
    cron.write()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculation of data for PDB.')
    parser.add_argument(
        '-l',
        '--log',
        help='Log file name.',
        default='updater.log')

    parser.add_argument(
        '-at',
        '--attemptPdb',
        help='Attempt Pdb coefficient. The number of tries for each pdb.',
        default='2')

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

    initial() # Remove the task relate to script_pdb
    schedule(
        args.log,
        args.attemptPdb,
        args.cathComplex, 
        args.cathChain, 
        args.cathDomain, 
        args.atomic, 
        args.initialAtomic, 
        args.amount, 
        args.atomicEmd, 
        args.hour)