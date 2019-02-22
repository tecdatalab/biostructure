from classes.Type_descriptor import Type_descriptor
import os

'''
Created on 22 feb. 2019

@author: luis98
'''
directory = "../db/"
directories_to_ignore = ["img","xml"]

def clean_directories(directoys):
    for k in directories_to_ignore:
        directoys.remove(k)
        
    i = 0
    while i < len(directoys):
        if(directoys[i].find('.') != -1):
            directoys.pop(i)
            i-=1
        else:
            i+=1
    return directoys

def name_for_directories(directoys,rootDir):
    result = []
    
    for directory in directoys:
        for dirName, _subdirList, fileList in os.walk(rootDir+directory+"/"):
            if fileList != []:
                add = dirName.replace(rootDir,"")
                add = add.replace("/","->")
                result.append(add)
                
    return result

def get_type_descriptors():
    result = []
    
    clean = clean_directories(list(os.listdir(directory)))
    names = name_for_directories(clean,directory)
    for i in names:
        add = Type_descriptor(0,i,i+" is a type of descriptor.")
        result.append(add)
    return result