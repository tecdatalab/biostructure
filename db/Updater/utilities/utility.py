import psutil, os
from csv import writer

'''
Created on 20 jun. 2020
@author: dnnxl

Last modified by @dnnxl
Date on 20 jun 2020
'''

# Return the memory (in MB) of a python process.
def memory():
    mem = psutil.Process(os.getpid()).memory_info().rss
    return (mem/1024**2)

# Append data in csv to a file.
def append_list_as_row(file_name, list_of_elem):
    with open(file_name, 'a+', newline='') as write_obj:
        csv_writer = writer(write_obj)
        csv_writer.writerow(list_of_elem)