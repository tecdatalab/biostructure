'''
Created on 28 abr. 2019

@author: luis98
'''
import os
import requests

dir = ""
emd_url = "http://ftp.ebi.ac.uk/pub/databases/emdb"

def download_file(emd_id):
    URL = "{1}/structures/EMD-{0}/map/emd_{0}.map.gz".format(emd_id,emd_url)
    response = requests.get(URL)
    with open('{0}temp/emd_{1}.map.gz'.format(dir,emd_id), 'wb') as file:
            file.write(response.content)
    os.system("gzip -d {0}temp/emd_{1}.map.gz".format(dir,emd_id))

def get_min_max_density(emd_id, contour):
    text = os.popen('./{0}em_volume {1} {2}temp/emd_{3}.map'.format(dir,contour,dir,emd_id)).read()
    list_a = text.split()
    return (float(list_a[10]),float(list_a[11]))

def generate_descriptor_file(emd_id, contour, exit_file_name):
    os.system("./{0}map2zernike {0}temp/emd_{1}.map -n 20 -c {2}  -p {0}temp/{3}".format(dir,emd_id,contour,exit_file_name))

def generate_descriptor_file_for_union(emd_id, contour_list, exit_file_name):
    temp_file_name = "temp"
    for i in contour_list:
        os.system("./{0}map2zernike {0}temp/emd_{1}.map -n 20 -c {2}  -p {0}temp/{3}".format(dir,emd_id,i,temp_file_name))
        os.system("tail -$(wc -l {0}temp/{1}.inv | awk '{{print $1-1}}') {0}temp/{1}.inv >> {0}temp/{2}.inv".format(dir,temp_file_name,exit_file_name))

    os.system("rm {0}temp/{1}.inv".format(dir,temp_file_name))

def generate_descriptors_files(emd_id, contour, std):
    max_density = get_min_max_density(emd_id,contour)[1]
    one_third_contour = contour + ((max_density-contour)/3)
    two_thirds_contour = contour + ((max_density-contour)*(2/3))
    one_std_contour = contour + std
    
    generate_descriptor_file(emd_id,contour,"contour")
    generate_descriptor_file_for_union(emd_id,[contour,one_third_contour],"one_third_contour")
    generate_descriptor_file_for_union(emd_id,[contour,two_thirds_contour],"two_thirds_contour")
    generate_descriptor_file_for_union(emd_id,[contour,one_third_contour,two_thirds_contour],"one_third_two_thirds_contour")
    generate_descriptor_file_for_union(emd_id,[contour,one_std_contour],"one_std_contour")

def remove_map(emd_id):
    os.system("rm {0}temp/emd_{1}.map".format(dir,emd_id))

def get_emd_descriptors(emd_id, contour, std):
    
    files = ["contour","one_third_contour","two_thirds_contour","one_third_two_thirds_contour","one_std_contour"]
    result = []
    
    generate_descriptors_files(emd_id, contour, std)
    
    for i in files:
        with open("{0}temp/{1}.inv".format(dir,i)) as f:
            array = [ float(line) for line in f]
            result.append(array)

    for i in files:
        os.system("rm {0}temp/{1}.inv".format(dir,i))
    
    return result



#download_file("0001")
#print(get_min_max_density("0001",0.018))
#generate_descriptor_file("0001",0.018,"normal_contour")
#result = get_emd_descriptors("0001",0.018, 0.0077065425)

#for i in result:
#    print(i)
