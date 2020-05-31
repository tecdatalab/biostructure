'''
Created on 10 nov. 2019
@author: luis98

Last modified: 21 may 2020
@By: dnnxl
'''
import os

class Calculate:
    
    def calculate_descriptors(self, atoms):
        dirpath = os.getcwd()
        foldername = os.path.basename(dirpath)
        if(foldername != "LZerD"):
            os.chdir("../generators/LZerD/")

        f= open("temp.pdb","w+")
        f.write(atoms)
        f.close()
        os.system("./runvdock_2.sh temp.pdb")
        #Get descriptors--------------
        f = open("temp_01.inv","r+")
        f.close()
        #---------------------------
        os.remove("temp.pdb")
        os.remove("temp_01.inv")
        os.chdir("../../updater")
        return []