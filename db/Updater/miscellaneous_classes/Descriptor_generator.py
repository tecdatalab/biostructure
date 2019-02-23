from classes.Descriptor import Descriptor
import os

'''
Created on 22 feb. 2019

@author: luis98
'''
directory = "../db/"

def get_descriptor(emd_id, typeDescriptor):
    result_temp = []
    
    path = directory + typeDescriptor.name
    path = path.replace("->","/")
    for i in os.listdir(path):
        if (i.find(str(emd_id)) != -1):
            file = open(path+ "/" +i, "r")
            for line in file.readlines():
                result_temp.append(line.replace('\n',''))
            file.close()
            
    return Descriptor(emd_id, None, result_temp, typeDescriptor.name)
    
def get_descriptors(emd_id, typeDescriptorList):
    result = []
    for descriptor in typeDescriptorList:
        result.append(get_descriptor(emd_id, descriptor))
    
    return result
    