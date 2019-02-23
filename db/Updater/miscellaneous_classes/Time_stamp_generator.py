from classes.Time_stamp import Time_stamp

'''
Created on 22 feb. 2019

@author: luis98
'''

directory = "../db/"

def get_time_stamps():
    result = []
    
    file = open(directory+"emdb_timestamps.txt", "r")
    for line in file.readlines()[1:]:
        line = line.replace("\n","")
        temp = line.split("\t")
        result.append(Time_stamp(temp[0],temp[1],temp[4],temp[5],temp[6]))
    file.close()
    return result
    
    