'''
Created on 18 mar. 2019

@author: luis98
'''
from classes.Emd_entry import Emd_entry
import pytest
from pytest_postgresql import factories  # Import the ability to use postgresql

emd_entry = Emd_entry(None, None, None, None, None, None,
                      None, None, None, None, None, None,
                      None, None, None, None, None, None,
                      None, None, None, None, None, None,
                      None, None, None, None, None, None,
                      None, None, None, None, None, None,
                      None, None, None, None, None, None,
                      None)


def test_0_set_get_id():

    id = 1
    emd_entry.id = id
    assert emd_entry.id == id


def test_1_set_get_full_name():

    full_name = "RsgA-30S ribosome-GMPPNP complex"
    emd_entry.full_name = full_name
    assert emd_entry.full_name == full_name


def test_2_set_get_acronym():

    acronym = "RsgA-30S"
    emd_entry.acronym = acronym
    assert emd_entry.acronym == acronym


def test_3_set_get_volume():

    volume = 0.1
    emd_entry.volume = volume
    assert emd_entry.volume == volume


def test_4_set_get_resolution():

    resolution = 0.1
    emd_entry.resolution = resolution
    assert emd_entry.resolution == resolution


def test_5_set_get_image_url():

    image_url = "imageURL"
    emd_entry.image_url = image_url
    assert emd_entry.image_url == image_url


def test_6_set_get_xml_url():

    xml_url = "xmlURL"
    emd_entry.xml_url = xml_url
    assert emd_entry.xml_url == xml_url


def test_7_set_get_map_url():

    map_url = "mapURL"
    emd_entry.map_url = map_url
    assert emd_entry.map_url == map_url


def test_8_set_get_map_id():

    map_id = 1
    emd_entry.map_id = map_id
    assert emd_entry.map_id == map_id


def test_9_set_get_map_file_information():

    file_information = {
        'format': 'CCP4',
        'sizeKb': 8584,
        'type': 'map',
        'file': 'emd_1364.map.gz'}
    emd_entry.file_information = file_information
    assert emd_entry.file_information == file_information


def test_10_set_get_map_data_type():

    data_type = "Image stored as Reals"
    emd_entry.data_type = data_type
    assert emd_entry.data_type == data_type


def test_11_set_get_map_num_columns():

    num_columns = 130
    emd_entry.num_columns = num_columns
    assert emd_entry.num_columns == num_columns


def test_12_set_get_map_num_rows():

    num_rows = 130
    emd_entry.num_rows = num_rows
    assert emd_entry.num_rows == num_rows


def test_13_set_get_map_num_sections():

    num_sections = 130
    emd_entry.num_sections = num_sections
    assert emd_entry.num_sections == num_sections


def test_14_set_get_map_origin_col():

    origin_col = -65
    emd_entry.origin_col = origin_col
    assert emd_entry.origin_col == origin_col


def test_15_set_get_map_origin_row():

    origin_row = -65
    emd_entry.origin_row = origin_row
    assert emd_entry.origin_row == origin_row


def test_16_set_get_map_origin_sec():

    origin_sec = -65
    emd_entry.origin_sec = origin_sec
    assert emd_entry.origin_sec == origin_sec


def test_17_set_get_map_limit_col():

    limit_col = 64
    emd_entry.limit_col = limit_col
    assert emd_entry.limit_col == limit_col


def test_18_set_get_map_limit_row():

    limit_row = 64
    emd_entry.limit_row = limit_row
    assert emd_entry.limit_row == limit_row


def test_19_set_get_map_limit_sec():

    limit_sec = 64
    emd_entry.limit_sec = limit_sec
    assert emd_entry.limit_sec == limit_sec


def test_20_set_get_map_spacing_col():

    spacing_col = 130
    emd_entry.spacing_col = spacing_col
    assert emd_entry.spacing_col == spacing_col


def test_21_set_get_map_spacing_row():

    spacing_row = 130
    emd_entry.spacing_row = spacing_row
    assert emd_entry.spacing_row == spacing_row


def test_22_set_get_map_spacing_sec():

    spacing_sec = 130
    emd_entry.spacing_sec = spacing_sec
    assert emd_entry.spacing_sec == spacing_sec


def test_23_set_get_map_cell_a():

    cell_a = {'units': 'A', 'value': 366.6}
    emd_entry.cell_a = cell_a
    assert emd_entry.cell_a == cell_a


def test_24_set_get_map_cell_b():

    cell_b = {'units': 'A', 'value': 366.6}
    emd_entry.cell_b = cell_b
    assert emd_entry.cell_b == cell_b


def test_25_set_get_map_cell_c():

    cell_c = {'units': 'A', 'value': 366.6}
    emd_entry.cell_c = cell_c
    assert emd_entry.cell_c == cell_c


def test_26_set_get_map_cell_alpha():

    cell_alpha = {'units': 'degrees', 'value': 90}
    emd_entry.cell_alpha = cell_alpha
    assert emd_entry.cell_alpha == cell_alpha


def test_27_set_get_map_cell_beta():

    cell_beta = {'units': 'degrees', 'value': 90}
    emd_entry.cell_beta = cell_beta
    assert emd_entry.cell_beta == cell_beta


def test_28_set_get_map_cell_gamma():

    cell_gamma = {'units': 'degrees', 'value': 90}
    emd_entry.cell_gamma = cell_gamma
    assert emd_entry.cell_gamma == cell_gamma


def test_29_set_get_map_axis_order_fast():

    axis_order_fast = "X"
    emd_entry.axis_order_fast = axis_order_fast
    assert emd_entry.axis_order_fast == axis_order_fast


def test_30_set_get_map_axis_order_medium():

    axis_order_medium = "Y"
    emd_entry.axis_order_medium = axis_order_medium
    assert emd_entry.axis_order_medium == axis_order_medium


def test_31_set_get_map_axis_order_slow():

    axis_order_slow = "Z"
    emd_entry.axis_order_slow = axis_order_slow
    assert emd_entry.axis_order_slow == axis_order_slow


def test_32_set_get_map_minimum():

    minimum = -78.5705
    emd_entry.minimum = minimum
    assert emd_entry.minimum == minimum


def test_33_set_get_map_maximum():

    maximum = 156.588
    emd_entry.maximum = maximum
    assert emd_entry.maximum == maximum


def test_34_set_get_map_average():

    average = -0.00988095
    emd_entry.average = average
    assert emd_entry.average == average


def test_35_set_get_map_std():

    std = 3.62791
    emd_entry.std = std
    assert emd_entry.std == std


def test_36_set_get_map_space_group_number():

    space_group_number = 1
    emd_entry.space_group_number = space_group_number
    assert emd_entry.space_group_number == space_group_number


def test_37_set_get_map_details():

    details = "::::EMDATABANK.org::::EMD-1364::::"
    emd_entry.details = details
    assert emd_entry.details == details


def test_38_set_get_map_pixel_x():

    pixel_x = {'units': 'A', 'value': 2.82}
    emd_entry.pixel_x = pixel_x
    assert emd_entry.pixel_x == pixel_x


def test_39_set_get_map_pixel_y():

    pixel_y = {'units': 'A', 'value': 2.82}
    emd_entry.pixel_y = pixel_y
    assert emd_entry.pixel_y == pixel_y


def test_40_set_get_map_pixel_z():

    pixel_z = {'units': 'A', 'value': 2.82}
    emd_entry.pixel_z = pixel_z
    assert emd_entry.pixel_z == pixel_z


def test_41_set_get_map_countour_level():

    countour_level = 9.07
    emd_entry.countour_level = countour_level
    assert emd_entry.countour_level == countour_level


def test_42_set_get_map_annotation_details():

    annotation_details = "Cryo-EM map of E.coli 70S ribosome"
    emd_entry.annotation_details = annotation_details
    assert emd_entry.annotation_details == annotation_details


def test_43_set_get_map_annotation_details():

    png_img_3d = "www.url.com"
    emd_entry.png_img_3d = png_img_3d
    assert emd_entry.png_img_3d == png_img_3d


def test_44_set_get_map_annotation_details():

    gif_img_3d = "www.url.com"
    emd_entry.gif_img_3d = gif_img_3d
    assert emd_entry.gif_img_3d == gif_img_3d


def test_45_insert_db(postgresql):

    create_map_information_table(postgresql)
    create_emd_entry_table(postgresql)
    emd_entrys = create_emd_entrys()
    insert_emd_entrys(postgresql, emd_entrys)
    result = get_emd_entrys(postgresql)
    assert result[0] == Emd_entry(
        1, "RsgA-30S ribosome-GMPPNP complex", "RsgA-30S",
        0.1, 0.1, "imageURL", "png_img_3d", "gif_img_3d", "xmlURL", "mapURL", 1,
        {'format': 'CCP4', 'sizeKb': 8584, 'type': 'map', 'file': 'emd_1364.map.gz'},
        "Image stored as Reals", 130, 130, 130, -65, -65, -65, 64, 64, 64, 130, 130, 130,
        {'units': 'A', 'value': 366.6}, {'units': 'A', 'value': 366.6},
        {'units': 'A', 'value': 366.6}, {'units': 'degrees', 'value': 90},
        {'units': 'degrees', 'value': 90}, {'units': 'degrees', 'value': 90},
        "X", "Y", "Z", -78.5705, 156.588, -0.00988095, 3.62791, 1,
        "::::EMDATABANK.org::::EMD-1364::::", {'units': 'A', 'value': 2.82},
        {'units': 'A', 'value': 2.82}, {'units': 'A', 'value': 2.82}, 9.07,
        "Cryo-EM map of E.coli 70S ribosome")


def test_46_update_db(postgresql):

    create_map_information_table(postgresql)
    create_emd_entry_table(postgresql)
    emd_entrys = create_emd_entrys()
    insert_emd_entrys(postgresql, emd_entrys)
    emd_entrys[0] = Emd_entry(
        1, "RsgA-30S ribosome-GMPPNP complexUp", "RsgA-30SUp",
        0.2, 0.2, "imageURLUp", "png_img_3dUp", "gif_img_3dUp", "xmlURLUp", "mapURLUp", 1,
        {'format': 'CCP4Up', 'sizeKb': 85842, 'type': 'mapUp', 'file': 'emd_1364.map.gzUp'},
        "Image stored as RealsUp", 230, 230, 230, -25, -25, -25, 24, 24, 24, 230, 230, 230,
        {'units': 'B', 'value': 366.2}, {'units': 'B', 'value': 366.2},
        {'units': 'B', 'value': 366.2}, {'units': 'degreesUP', 'value': 92},
        {'units': 'degreesUp', 'value': 92}, {'units': 'degreesUp', 'value': 92},
        "Z", "X", "Y", -78.5702, 156.582, -0.00988092, 3.62792, 2,
        "::::EMDATABANK.org::::EMD-1364::::Up", {'units': 'B', 'value': 2.83},
        {'units': 'B', 'value': 2.32}, {'units': 'B', 'value': 2.32}, 9.37,
        "Cryo-EM map of E.coli 70S ribosomeUp")

    update_emd_entrys(postgresql, emd_entrys)
    result = get_emd_entrys(postgresql)
    assert result[0] == Emd_entry(
        1, "RsgA-30S ribosome-GMPPNP complexUp", "RsgA-30SUp",
        0.2, 0.2, "imageURLUp", "png_img_3dUp", "gif_img_3dUp", "xmlURLUp", "mapURLUp", 1,
        {'format': 'CCP4Up', 'sizeKb': 85842, 'type': 'mapUp', 'file': 'emd_1364.map.gzUp'},
        "Image stored as RealsUp", 230, 230, 230, -25, -25, -25, 24, 24, 24, 230, 230, 230,
        {'units': 'B', 'value': 366.2}, {'units': 'B', 'value': 366.2},
        {'units': 'B', 'value': 366.2}, {'units': 'degreesUP', 'value': 92},
        {'units': 'degreesUp', 'value': 92}, {'units': 'degreesUp', 'value': 92},
        "Z", "X", "Y", -78.5702, 156.582, -0.00988092, 3.62792, 2,
        "::::EMDATABANK.org::::EMD-1364::::Up", {'units': 'B', 'value': 2.83},
        {'units': 'B', 'value': 2.32}, {'units': 'B', 'value': 2.32}, 9.37,
        "Cryo-EM map of E.coli 70S ribosomeUp")


def get_emd_entrys(connection):

    cursor = connection.cursor()
    cursor.execute(
        "SELECT * FROM emd_entry em INNER JOIN map_information map ON (em.map_information_id = map.id)")
    emd_entrys = [
        Emd_entry(
            record[0],
            record[1],
            record[2],
            record[3],
            record[4],
            record[5],
            record[6],
            record[7],
            record[8],
            record[9],
            record[10],
            record[12],
            record[13],
            record[14],
            record[15],
            record[16],
            record[17],
            record[18],
            record[19],
            record[20],
            record[21],
            record[22],
            record[23],
            record[24],
            record[25],
            record[26],
            record[27],
            record[28],
            record[29],
            record[30],
            record[31],
            record[32],
            record[33],
            record[34],
            record[35],
            record[36],
            record[37],
            record[38],
            record[39],
            record[40],
            record[41],
            record[42],
            record[43],
            record[44],
            record[45]) for record in cursor]
    cursor.close()

    return emd_entrys


def insert_emd_entrys(connection, emd_entrys):
    cursor = connection.cursor()

    for emd_entry in emd_entrys:
        emd_entry.insert_db(cursor)

    connection.commit()
    cursor.close()


def update_emd_entrys(connection, emd_entrys):
    cursor = connection.cursor()

    for emd_entry in emd_entrys:
        emd_entry.update_db(cursor)

    connection.commit()
    cursor.close()


def create_emd_entrys():
    resultado = []
    resultado.append(Emd_entry(
        1, "RsgA-30S ribosome-GMPPNP complex", "RsgA-30S",
        0.1, 0.1, "imageURL", "png_img_3d", "gif_img_3d", "xmlURL", "mapURL", 1,
        {'format': 'CCP4', 'sizeKb': 8584, 'type': 'map', 'file': 'emd_1364.map.gz'},
        "Image stored as Reals", 130, 130, 130, -65, -65, -65, 64, 64, 64, 130, 130, 130,
        {'units': 'A', 'value': 366.6}, {'units': 'A', 'value': 366.6},
        {'units': 'A', 'value': 366.6}, {'units': 'degrees', 'value': 90},
        {'units': 'degrees', 'value': 90}, {'units': 'degrees', 'value': 90},
        "X", "Y", "Z", -78.5705, 156.588, -0.00988095, 3.62791, 1,
        "::::EMDATABANK.org::::EMD-1364::::", {'units': 'A', 'value': 2.82},
        {'units': 'A', 'value': 2.82}, {'units': 'A', 'value': 2.82}, 9.07,
        "Cryo-EM map of E.coli 70S ribosome"))
    return resultado


def create_emd_entry_table(connection):

    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE emd_entry(\
        id INT PRIMARY KEY,\
        full_name TEXT NOT NULL,\
        acronym TEXT NOT NULL,\
        volume FLOAT8 NOT NULL,\
        resolution FLOAT8 NOT NULL,\
        image_url TEXT,\
        png_img_3d TEXT,\
        gif_img_3d TEXT,\
        xml_url TEXT NOT NULL,\
        map_url TEXT NOT NULL,\
        map_information_id INT REFERENCES map_information(id) NOT NULL\
        );")
    connection.commit()
    cursor.close()


def create_map_information_table(connection):

    cursor = connection.cursor()
    cursor.execute(
        "CREATE TABLE map_information(\
        id SERIAL PRIMARY KEY,\
        file_information JSON NOT NULL,\
        data_type TEXT NOT NULL,\
        num_columns INT NOT NULL,\
        num_rows INT NOT NULL,\
        num_sections INT NOT NULL,\
        origin_col INT NOT NULL,\
        origin_row INT NOT NULL,\
        origin_sec INT NOT NULL,\
        limit_col INT NOT NULL,\
        limit_row INT NOT NULL,\
        limit_sec INT NOT NULL,\
        spacing_col INT NOT NULL,\
        spacing_row INT NOT NULL,\
        spacing_sec INT NOT NULL,\
        cell_a JSON NOT NULL,\
        cell_b JSON NOT NULL,\
        cell_c JSON NOT NULL,\
        cell_alpha JSON NOT NULL,\
        cell_beta JSON NOT NULL,\
        cell_gamma JSON NOT NULL,\
        axis_order_fast CHAR(1) NOT NULL,\
        axis_order_medium CHAR(1) NOT NULL,\
        axis_order_slow CHAR(1) NOT NULL,\
        minimum FLOAT8 NOT NULL,\
        maximum FLOAT8 NOT NULL,\
        average FLOAT8 NOT NULL,\
        std FLOAT8 NOT NULL,\
        space_group_number INT NOT NULL,\
        details TEXT NOT NULL,\
        pixel_x JSON NOT NULL,\
        pixel_y JSON NOT NULL,\
        pixel_z JSON NOT NULL,\
        countour_level FLOAT8 NOT NULL,\
        annotation_details TEXT\
        );")
    connection.commit()
    cursor.close()
