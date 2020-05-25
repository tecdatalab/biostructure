'''
file logging.py
@author Danny Xie 
This file will keep all the functions to log a file
'''

import logging, sys, os
from constants.constants import *

global logging

class Log(object):

    __file = None

    def __init__(self, fileName):
        self.__file = fileName

    # Function init_log_file 
    # Description: The following function init a log file. with the name of the logFileName, the basic configuration.
    # By: Danny Xie Li
    # Created: 16/05/2020
    # Last modification: 16/05/2020

    def init_log_file(self, file):
        Constant.config['file'] = file
        logging.basicConfig(
            level=logging.INFO,
            format = Constant.formatLog,
            handlers=[
                logging.FileHandler(self.__file),
                logging.StreamHandler()
            ]
        )

    # Function generate_info_message 
    # Description: The following function generate an info in the log file and in the console.
    # By: Danny Xie Li
    # Created: 16/05/2020
    # Last modification: 16/05/2020

    def generate_info_message(self, personalized_info, compute_info):
        Constant.config['type'] = str(compute_info)
        Constant.config['error'] = ''
        logging.info(personalized_info, extra = Constant.config)

    # Function generate_error_message 
    # Description: The following function generate an error in the log file and in the console.
    # By: Danny Xie Li
    # Created: 16/05/2020
    # Last modification: 16/05/2020

    def generate_error_message(self, personalized_error, compute_error, error):
        Constant.config['type'] = str(compute_error)
        Constant.config['error'] = str(error)
        logging.error(personalized_error, extra = Constant.config)
