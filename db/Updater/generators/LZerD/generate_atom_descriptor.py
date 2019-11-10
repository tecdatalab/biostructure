'''
Created on 10 nov. 2019

@author: luis98
'''
import os

class Calculate:
    
    def calculate_descriptors(self, atoms):
        os.chdir("../generators/LZerD")
        f= open("temp.pdb","w+")
        f.write(atoms)
        f.close()
        os.system("./runvdock_2.sh temp.pdb");
        #Get descriptors--------------
        f = open("temp_01.inv","r+")
        f.close()
        #---------------------------
        os.remove("temp.pdb")
        os.remove("temp_01.inv")
        os.chdir("../../updater")
        return []