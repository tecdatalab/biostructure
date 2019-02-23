from classes.Update import Update

'''
Created on 22 feb. 2019

@author: luis98
'''

directory = "../db/"

def get_updates():
    result = []
    
    file = open(directory+"latest.updated.date", "r")
    for line in file.readlines():
        line = line.replace("\n","")
        result.append(Update(line))
    file.close()
    return result