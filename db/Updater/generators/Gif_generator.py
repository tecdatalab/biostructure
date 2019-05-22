'''
Created on 22 may. 2019

@author: luis98
'''
from PIL import Image
import glob, os

dir = ""

def generateGif(emd_id):
    frames = []
    
    for file in os.listdir("{0}temp_images".format(dir)):
        image_dir = dir + "temp_images/" + file
        frames.append(Image.open(image_dir))
    frames[0].save('{0}db_images/{1}.gif'.format(dir, emd_id), format='GIF', append_images=frames[1:], save_all=True, duration=100, loop=0)
    

#generateGif(12)