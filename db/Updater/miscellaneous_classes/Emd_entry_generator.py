from classes.Emd_entry import Emd_entry

'''
Created on 22 feb. 2019

@author: luis98
'''

directory = "../db/"

def get_second_element_file(fileName,emd_id):
    file = open(directory+fileName, "r")
    for line in file.readlines():
        if line.split("\t")[0] == str(emd_id):
            file.close()
            return line.split("\t")[1].replace("\n","")
    file.close()

def get_full_name(emd_id):
    return get_second_element_file("fullname.txt",emd_id)
    
def get_acronym(emd_id):
    return get_second_element_file("name.txt",emd_id)
    
def get_volume(emd_id):
    return get_second_element_file("volume.txt",emd_id)
    
def get_resolution(emd_id):
    return get_second_element_file("resolutions.txt",emd_id)

def get_emd_entry(emd_id):
    emd_id = str(emd_id)
    full_name = get_full_name(emd_id)
    acronym =  get_acronym(emd_id)
    volume =  get_volume(emd_id)
    resolution = get_resolution(emd_id)
    
    return Emd_entry(emd_id,full_name,acronym,volume,resolution,
                     None,None,None,None)
def get_emd_entries():
    result = []
    
    file = open(directory+"idlist.txt", "r")
    for line in file.readlines():
        result.append((line.replace("\n","")))
    file.close()
    return result
    
    
    