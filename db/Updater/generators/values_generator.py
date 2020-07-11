'''
Created on 28 abr. 2019

@author: luis98

Last modification on 14 jun. 2020
@author: dnnxl
'''
import os
import os.path
import requests
from clint.textui import progress

dir = ""

def download_file(emd_id):
    URL = "http://ftp.wwpdb.org/pub/emdb/structures/EMD-{0}/map/emd_{0}.map.gz".format(emd_id)
    response = requests.get(URL, stream=True)
    with open('{0}temp/emd_{1}.map.gz'.format(dir, emd_id), 'wb') as file:
        total_length = int(response.headers.get('content-length'))
        for chunk in progress.bar(response.iter_content(chunk_size = 8192), expected_size=(total_length/8192) + 1):
            if chunk:
                file.write(chunk)
                file.flush()
                os.fsync(file.fileno())
    os.system("gzip -d -f {0}temp/emd_{1}.map.gz".format(dir, emd_id))


def get_min_max_density(emd_id, contour):
    text = os.popen(
        './{0}em_volume {1} {2}temp/emd_{3}.map'.format(dir, contour, dir, emd_id)).read()
    list_a = text.split()
    return (float(list_a[10]), float(list_a[11]))


def generate_descriptor_file(emd_id, contour, exit_file_name):
    execute_string = "./{0}map2zernike {0}temp/emd_{1}.map -n 20 -c {2} -p {0}temp/{3} ".format(
            dir,
            emd_id,
            contour,
            exit_file_name)
    os.system(
        execute_string)
    if (not os.path.exists("{0}temp/{1}.inv".format(dir, exit_file_name))):
        raise Exception('Emd descriptor file was not created.')


def generate_descriptor_file_for_union(emd_id, contour_list, exit_file_name):
    temp_file_name = "temp"
    for i in contour_list:
        generate_descriptor_file(emd_id, i, temp_file_name)
        os.system(
            "tail -$(wc -l {0}temp/{1}.inv | awk '{{print $1-1}}') {0}temp/{1}.inv >> {0}temp/{2}.inv ".format(
                dir,
                temp_file_name,
                exit_file_name))

    os.system("rm {0}temp/{1}.inv ".format(dir, temp_file_name))


def generate_descriptors_files(emd_id, contour, std):
    max_density = get_min_max_density(emd_id, contour)[1]
    one_third_contour = contour + ((max_density - contour) / 3)
    two_thirds_contour = contour + ((max_density - contour) * (2 / 3))
    one_std_contour = contour + std

    generate_descriptor_file(emd_id, contour, "contour")
    generate_descriptor_file_for_union(
        emd_id, [contour, one_third_contour], "one_third_contour")
    generate_descriptor_file_for_union(
        emd_id, [contour, two_thirds_contour], "two_thirds_contour")
    generate_descriptor_file_for_union(
        emd_id, [
            contour, one_third_contour, two_thirds_contour], "one_third_two_thirds_contour")
    generate_descriptor_file_for_union(
        emd_id, [contour, one_std_contour], "one_std_contour")


def remove_map(emd_id):
    os.system("rm -rf {0}temp ".format(dir))
    os.system("mkdir {0}temp ".format(dir))


def get_emd_descriptors(emd_id, contour, std):

    files = [
        "contour",
        "one_third_contour",
        "two_thirds_contour",
        "one_third_two_thirds_contour",
        "one_std_contour"]
    result = []

    generate_descriptors_files(emd_id, contour, std)

    for i in files:
        with open("{0}temp/{1}.inv".format(dir, i)) as f:
            array = [float(line) for line in f]
            result.append(array)

    for i in files:
        os.system("rm {0}temp/{1}.inv ".format(dir, i))

    return result


# download_file("0001")
# print(get_min_max_density("0001",0.018))
# generate_descriptor_file("0001",0.018,"normal_contour")
#result = get_emd_descriptors("0001",0.018, 0.0077065425)

# for i in result:
#    print(i)
