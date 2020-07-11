import sys, argparse, os
sys.path.append('../')
import datetime

from connections.sql_connection import SQL_connection
from classes.update import Update
from classes.binnacle import Binnacle
from main_emd import main

def main_script_emd(connec_sql, mode, image, descriptor, startEmd, amount):
    cursor_sql = connec_sql.get_cursor()

    update_temp = Update(datetime.datetime.now())
    update_temp.insert_update_db(connec_sql, cursor_sql)
    connec_sql.commit()

    temp_binnacle = Binnacle()
    temp_binnacle.set_last_update(cursor_sql)
    temp_emd_id = temp_binnacle.get_last_emd_id_update(cursor_sql)
    
    connec_sql.close_connection()

    if(temp_emd_id == None):
        main(startEmd, mode, image, descriptor, startEmd + amount)
    else:
        main(temp_emd_id, mode, image, descriptor, temp_emd_id + amount)

def init_script_emd(mode, image, descriptor, startEmd, amount):
    connec_sql = SQL_connection()
    connec_sql.init_connection()
    main_script_emd(connec_sql, mode, image, descriptor, int(startEmd), int(amount))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Scheduling the execution of updater')
    parser.add_argument(
        '-j',
        '--json',
        help='Json file name.',
        default='scheduling.json')
    parser.add_argument(
        '-ie',
        '--initialEMD',
        help='Initial EMD for execution.',
        default='0001')

    parser.add_argument(
        '-a',
        '--amount',
        help='Amount of emd to execute.',

        default='1')
    
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
    filename = args.json
    init_script_emd(args.mode, args.image, args.descriptor, args.initialEMD, args.amount)