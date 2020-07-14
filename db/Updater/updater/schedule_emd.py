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

def schedule(log, mode, image, descriptor, startEmd, hour, amount):
    cron = CronTab(user=True)
    temp_command = "sudo python3 {}/script_emd.py -l {} -m {} -i {} -d {} -ie {} -a {}".format(directory, log, mode, image, descriptor, startEmd, amount)
    print(temp_command)
    job = cron.new(command=temp_command)
    job.every(hour).hours()
    cron.write()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Scheduling the execution of updater')

    parser.add_argument(
        '-l',
        '--log',
        help='Log file name.',
        default='updater.log')
    parser.add_argument(
        '-ho',
        '--hour',
        help='Time schedule.',
        default='1')
    parser.add_argument(
        '-a',
        '--amount',
        help='Amount of emd to execute.',
        default='1')

    parser.add_argument(
        '-ie',
        '--initialEMD',
        help='Initial EMD for execution.',
        default='0001')

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
    print(args)
    initial()
    schedule(args.log, args.mode, args.image, args.descriptor, args.initialEMD, args.hour, args.amount)