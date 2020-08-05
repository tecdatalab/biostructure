import psutil, os
from csv import writer

'''
Created on 20 jun. 2020
@author: dnnxl

Last modified by @dnnxl
Date on 5 Aug. 2020
'''

# Return the memory (in MB) of a python process.
def memory():
    '''
    Inputs:
        None.
    Output:
        The memory (in MB) of a python process.
    '''
    mem = psutil.Process(os.getpid()).memory_info().rss
    return (mem/1024**2)

# Append data in csv to a file.
def append_list_as_row(file_name, list_of_elem):
    '''
    Inputs:
        file_name     -  String, the name of the file name.
        list_of_elem  -  List, the list to append in the file.
    Output:
        None.
    '''
    with open(file_name, 'a+', newline='') as write_obj:
        csv_writer = writer(write_obj)
        csv_writer.writerow(list_of_elem)

# Format result of the get_volume_map to list of strings
def format_volume_map(volume_map):
    '''
    Inputs:
        volume_map  -  Dictionary of the volumes of the map with 
                       the following format 
                       {contour:volume, ... , contour_n:volume_n}
    Output:
        Return a list of strings with the recommended contour and 
        the volume values of the map with the following format:  
        ["{ contour : value_contour}", ... , "{volume_n : value_volume_n}"].
    '''
    formated = []
    for contour, volume in volume_map.items():
        elem_formatted = {"contour":contour, "volume":volume}
        formated = formated + [str(elem_formatted)]
    return formated