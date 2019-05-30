'''
Created on 22 may. 2019

@author: luis98
'''
from PIL import Image
import glob, os
import shutil


dir = ""
images_path = "assets/img/front/{0}.png"
gif_path = "assets/img/gif/{0}.gif"

def generateGif(emd_id):
    frames = []
    
    for file in os.listdir("{0}temp_images".format(dir)):
        image_dir = dir + "temp_images/" + file
        frames.append(Image.open(image_dir))
    frames[0].save('{0}gif/{1}.gif'.format(dir, emd_id), format='GIF', append_images=frames[1:], save_all=True, duration=100, loop=0)
    shutil.copy(dir + "temp_images/1.png", '{0}front/{1}.png'.format(dir, emd_id))
    
    
    return (images_path.format(emd_id),gif_path.format(emd_id))

generateGif(12)