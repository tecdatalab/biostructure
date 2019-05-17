'''
Created on 28 abr. 2019

@author: luis98
'''
import pytest
import os
from descriptors_generator import values_generator as valg
import filecmp


valg.dir = "../descriptors_generator/"
os.system("cp temp/emd_0001.map {0}/temp".format(valg.dir))

def test_0_get_min_max_density():
    result = valg.get_min_max_density("0001",0.018)
    assert result == (-0.107885, 0.226883)

def test_1_generate_descriptor_file():
    valg.generate_descriptor_file("0001",0.018,"contour")
    temp = filecmp.cmp("temp/contour.inv", "{0}temp/contour.inv".format(valg.dir))
    assert temp

def test_2_generate_descriptor_file_for_union():
    max_density = valg.get_min_max_density("0001",0.018)[1]
    one_third_contour = 0.018 + ((max_density-0.018)/3)
    valg.generate_descriptor_file_for_union("0001",[0.018,one_third_contour],"one_third_contour")
    temp = filecmp.cmp("temp/one_third_contour.inv", "{0}temp/one_third_contour.inv".format(valg.dir))
    assert temp
    
def test_3_generate_descriptors_files():
    os.system("rm -r {0}/temp/".format(valg.dir))
    os.system("mkdir {0}/temp".format(valg.dir))
    os.system("cp temp/emd_0001.map {0}/temp".format(valg.dir))    
    valg.generate_descriptors_files("0001",0.018, 0.0077065425)
    temp_0 = filecmp.cmp("temp/contour.inv", "{0}temp/contour.inv".format(valg.dir))
    temp_1 = filecmp.cmp("temp/one_third_contour.inv", "{0}temp/one_third_contour.inv".format(valg.dir))
    temp_2 = filecmp.cmp("temp/two_thirds_contour.inv", "{0}temp/two_thirds_contour.inv".format(valg.dir))
    temp_3 = filecmp.cmp("temp/one_third_two_thirds_contour.inv", "{0}temp/one_third_two_thirds_contour.inv".format(valg.dir))
    temp_4 = filecmp.cmp("temp/one_std_contour.inv", "{0}temp/one_std_contour.inv".format(valg.dir))
    assert temp_0 and temp_1 and temp_2 and temp_3 and temp_4

def test_n_clean():
    os.system("rm -r {0}/temp/".format(valg.dir))
    os.system("mkdir {0}/temp".format(valg.dir))
    assert True
