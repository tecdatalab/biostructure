'''
Created on 22 may. 2019
@author: luis98

Last modification on 6 aug. 2020
@author: dnnxl
'''
import sys
sys.path.append('../')
from PIL import Image
import glob
import os
import shutil
from generators.reader import Reader
from generators.visualizer import Visualizer
from glumpy import app
import os
from numpy import sort

dir = ""
images_path = "https://biostructure.s3.us-east-2.amazonaws.com/front/{0}.png"
gif_path = "https://biostructure.s3.us-east-2.amazonaws.com/gif/{0}.gif"


def generateGif(emd_id, levelp):
    mapReader = Reader()
    mapReader.open("{0}temp/emd_{1}.map".format(dir, emd_id))
    # Get map object
    myMap = mapReader.read()
    # Create visualizer with a map surface threshold level
    # Otherwise use otsu threshold
    v = Visualizer(myMap, level=levelp)
    #v = Visualizer(myMap)
    # Watershed
    v.segmentate()
    # add corresponding atomic structure
    # v.add_structure("pdb6gh5.ent")
    v.show(export=True, time=3, export_path=dir)

    frames = []

    dirs = os.listdir("{0}export".format(dir))
    dirs = sorted(dirs, key=lambda dir: int(dir.split(".")[0]))
    for file in dirs:
        image_dir = dir + "export/" + file
        frames.append(Image.open(image_dir))
    frames[0].save('{0}gif/{1}.gif'.format(dir, emd_id), format='GIF',
                   append_images=frames[1:], save_all=True, duration=166, loop=0)
    shutil.copy("{0}export/{1}".format(dir,
                                       dirs[0]), '{0}front/{1}.png'.format(dir, emd_id))

    shutil.rmtree("{0}export".format(dir))
    os.mkdir("{0}export".format(dir))
    del(v)
    del(myMap)
    del(mapReader)
    return (images_path.format(emd_id), gif_path.format(emd_id))

# generateGif(12)
