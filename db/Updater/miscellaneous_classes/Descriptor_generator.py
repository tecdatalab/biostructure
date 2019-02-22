from classes.Descriptor import Descriptor
import os

'''
Created on 22 feb. 2019

@author: luis98
'''
directory = "../db/"

def get_descriptor(emd_id, descriptor):
    result_temp = []
    
    path = directory + descriptor.name
    path = path.replace("->","/")
    for i in os.listdir(path):
        if (i.find(str(emd_id)) != -1):
            file = open(path+ "/" +i, "r")
            for line in file.readlines():
                result_temp.append(line.replace('\n',''))
            file.close()
            
    return Descriptor(emd_id, None, result_temp, descriptor.name)
    
def get_descriptors(emd_id, descriptorList):
    result = []
    for descriptor in descriptorList:
        result.append(get_descriptor(emd_id, descriptor))
    
    return result
    