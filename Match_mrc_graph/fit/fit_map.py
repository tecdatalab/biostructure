'''
Created on Jul 19, 2020

@author: luis98
'''
import os
from subprocess import check_output, CalledProcessError
from tempfile import TemporaryFile
import shutil
from fit.fitMapResult import FitMapResult

def get_out(*args):
    with TemporaryFile() as t:
        try:
            out = check_output(args, stderr=t)
            return  0, out
        except CalledProcessError as e:
            t.seek(0)
            return e.returncode, t.read()

#Fit map0 into map1
def fit_map_in_map(map0_path, map1_path, map0_level = None, map1_level = None):
    path = "./temp_map"
    try:
        shutil.rmtree(path)
    except:
        pass
    os.mkdir(path)
    f = open(path+"/fit.cmd","w+")
    f.write("fitmap #0 #1 \r\n")
    
    if map0_level !=None:
        f.write("volume #0 level "+str(map0_level).replace(".", ",")+"\r\n")
    
    if map1_level !=None:
        f.write("volume #1 level "+str(map1_level).replace(".", ",")+"\r\n")
    
    f.close()
    
    map0_real_path = os.path.abspath(map0_path)
    map1_real_path = os.path.abspath(map1_path)
    commands_real_path = os.path.abspath(path+"/fit.cmd")
    
    _error, exit = get_out("chimera","--nogui", map0_real_path, map1_real_path, commands_real_path)

    shutil.rmtree(path)
    text = exit.decode("utf-8") 
    return FitMapResult(text)
    
    
    