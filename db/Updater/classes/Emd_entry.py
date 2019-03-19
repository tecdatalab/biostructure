'''
Created on 22 feb. 2019

@author: luis98
'''
from psycopg2 import sql
import json

class Map_information(object):
    '''
    classdocs
    '''
    __map_id = None
    __map_file_information = None
    __map_data_type = None
    __map_num_columns = None
    __map_num_rows = None
    __map_num_sections = None
    __map_origin_col = None
    __map_origin_row = None
    __map_origin_sec = None
    __map_limit_col = None
    __map_limit_row = None
    __map_limit_sec = None
    __map_spacing_col = None
    __map_spacing_row = None
    __map_spacing_sec = None
    __map_cell_a = None
    __map_cell_b = None
    __map_cell_c = None
    __map_cell_alpha = None
    __map_cell_beta = None
    __map_cell_gamma = None
    __map_axis_order_fast = None
    __map_axis_order_medium = None
    __map_axis_order_slow = None
    __map_minimum = None
    __map_maximum = None
    __map_average = None
    __map_std = None
    __map_space_group_number = None
    __map_details = None
    __map_pixel_x = None
    __map_pixel_y = None
    __map_pixel_z = None
    __map_countour_level = None
    __map_annotation_details = None
    
    
    def __init__(self, map_id, map_file_information, map_data_type, map_num_columns, map_num_rows, map_num_sections, map_origin_col, map_origin_row, map_origin_sec, map_limit_col, map_limit_row, map_limit_sec, map_spacing_col, map_spacing_row, map_spacing_sec, map_cell_a, map_cell_b, map_cell_c, map_cell_alpha, map_cell_beta, map_cell_gamma, map_axis_order_fast, map_axis_order_medium, map_axis_order_slow, map_minimum, map_maximum, map_average, map_std, map_space_group_number, map_details, map_pixel_x, map_pixel_y, map_pixel_z, map_countour_level, map_annotation_details):
        self.__map_id = map_id
        self.__map_file_information = map_file_information
        self.__map_data_type = map_data_type
        self.__map_num_columns = map_num_columns
        self.__map_num_rows = map_num_rows
        self.__map_num_sections = map_num_sections
        self.__map_origin_col = map_origin_col
        self.__map_origin_row = map_origin_row
        self.__map_origin_sec = map_origin_sec
        self.__map_limit_col = map_limit_col
        self.__map_limit_row = map_limit_row
        self.__map_limit_sec = map_limit_sec
        self.__map_spacing_col = map_spacing_col
        self.__map_spacing_row = map_spacing_row
        self.__map_spacing_sec = map_spacing_sec
        self.__map_cell_a = map_cell_a
        self.__map_cell_b = map_cell_b
        self.__map_cell_c = map_cell_c
        self.__map_cell_alpha = map_cell_alpha
        self.__map_cell_beta = map_cell_beta
        self.__map_cell_gamma = map_cell_gamma
        self.__map_axis_order_fast = map_axis_order_fast
        self.__map_axis_order_medium = map_axis_order_medium
        self.__map_axis_order_slow = map_axis_order_slow
        self.__map_minimum = map_minimum
        self.__map_maximum = map_maximum
        self.__map_average = map_average
        self.__map_std = map_std
        self.__map_space_group_number = map_space_group_number
        self.__map_details = map_details
        self.__map_pixel_x = map_pixel_x
        self.__map_pixel_y = map_pixel_y
        self.__map_pixel_z = map_pixel_z
        self.__map_countour_level = map_countour_level
        self.__map_annotation_details = map_annotation_details


    def get_map_id(self):
        return self.__map_id


    def get_map_file_information(self):
        return self.__map_file_information


    def get_map_data_type(self):
        return self.__map_data_type


    def get_map_num_columns(self):
        return self.__map_num_columns


    def get_map_num_rows(self):
        return self.__map_num_rows


    def get_map_num_sections(self):
        return self.__map_num_sections


    def get_map_origin_col(self):
        return self.__map_origin_col


    def get_map_origin_row(self):
        return self.__map_origin_row


    def get_map_origin_sec(self):
        return self.__map_origin_sec


    def get_map_limit_col(self):
        return self.__map_limit_col


    def get_map_limit_row(self):
        return self.__map_limit_row


    def get_map_limit_sec(self):
        return self.__map_limit_sec


    def get_map_spacing_col(self):
        return self.__map_spacing_col


    def get_map_spacing_row(self):
        return self.__map_spacing_row


    def get_map_spacing_sec(self):
        return self.__map_spacing_sec


    def get_map_cell_a(self):
        return self.__map_cell_a


    def get_map_cell_b(self):
        return self.__map_cell_b


    def get_map_cell_c(self):
        return self.__map_cell_c


    def get_map_cell_alpha(self):
        return self.__map_cell_alpha


    def get_map_cell_beta(self):
        return self.__map_cell_beta


    def get_map_cell_gamma(self):
        return self.__map_cell_gamma


    def get_map_axis_order_fast(self):
        return self.__map_axis_order_fast


    def get_map_axis_order_medium(self):
        return self.__map_axis_order_medium


    def get_map_axis_order_slow(self):
        return self.__map_axis_order_slow


    def get_map_minimum(self):
        return self.__map_minimum


    def get_map_maximum(self):
        return self.__map_maximum


    def get_map_average(self):
        return self.__map_average


    def get_map_std(self):
        return self.__map_std


    def get_map_space_group_number(self):
        return self.__map_space_group_number


    def get_map_details(self):
        return self.__map_details


    def get_map_pixel_x(self):
        return self.__map_pixel_x


    def get_map_pixel_y(self):
        return self.__map_pixel_y


    def get_map_pixel_z(self):
        return self.__map_pixel_z


    def get_map_countour_level(self):
        return self.__map_countour_level


    def get_map_annotation_details(self):
        return self.__map_annotation_details


    def set_map_id(self, value):
        self.__map_id = value


    def set_map_file_information(self, value):
        self.__map_file_information = value


    def set_map_data_type(self, value):
        self.__map_data_type = value


    def set_map_num_columns(self, value):
        self.__map_num_columns = value


    def set_map_num_rows(self, value):
        self.__map_num_rows = value


    def set_map_num_sections(self, value):
        self.__map_num_sections = value


    def set_map_origin_col(self, value):
        self.__map_origin_col = value


    def set_map_origin_row(self, value):
        self.__map_origin_row = value


    def set_map_origin_sec(self, value):
        self.__map_origin_sec = value


    def set_map_limit_col(self, value):
        self.__map_limit_col = value


    def set_map_limit_row(self, value):
        self.__map_limit_row = value


    def set_map_limit_sec(self, value):
        self.__map_limit_sec = value


    def set_map_spacing_col(self, value):
        self.__map_spacing_col = value


    def set_map_spacing_row(self, value):
        self.__map_spacing_row = value


    def set_map_spacing_sec(self, value):
        self.__map_spacing_sec = value


    def set_map_cell_a(self, value):
        self.__map_cell_a = value


    def set_map_cell_b(self, value):
        self.__map_cell_b = value


    def set_map_cell_c(self, value):
        self.__map_cell_c = value


    def set_map_cell_alpha(self, value):
        self.__map_cell_alpha = value


    def set_map_cell_beta(self, value):
        self.__map_cell_beta = value


    def set_map_cell_gamma(self, value):
        self.__map_cell_gamma = value


    def set_map_axis_order_fast(self, value):
        self.__map_axis_order_fast = value


    def set_map_axis_order_medium(self, value):
        self.__map_axis_order_medium = value


    def set_map_axis_order_slow(self, value):
        self.__map_axis_order_slow = value


    def set_map_minimum(self, value):
        self.__map_minimum = value


    def set_map_maximum(self, value):
        self.__map_maximum = value


    def set_map_average(self, value):
        self.__map_average = value


    def set_map_std(self, value):
        self.__map_std = value


    def set_map_space_group_number(self, value):
        self.__map_space_group_number = value


    def set_map_details(self, value):
        self.__map_details = value


    def set_map_pixel_x(self, value):
        self.__map_pixel_x = value


    def set_map_pixel_y(self, value):
        self.__map_pixel_y = value


    def set_map_pixel_z(self, value):
        self.__map_pixel_z = value


    def set_map_countour_level(self, value):
        self.__map_countour_level = value


    def set_map_annotation_details(self, value):
        self.__map_annotation_details = value


    def del_map_id(self):
        del self.__map_id


    def del_map_file_information(self):
        del self.__map_file_information


    def del_map_data_type(self):
        del self.__map_data_type


    def del_map_num_columns(self):
        del self.__map_num_columns


    def del_map_num_rows(self):
        del self.__map_num_rows


    def del_map_num_sections(self):
        del self.__map_num_sections


    def del_map_origin_col(self):
        del self.__map_origin_col


    def del_map_origin_row(self):
        del self.__map_origin_row


    def del_map_origin_sec(self):
        del self.__map_origin_sec


    def del_map_limit_col(self):
        del self.__map_limit_col


    def del_map_limit_row(self):
        del self.__map_limit_row


    def del_map_limit_sec(self):
        del self.__map_limit_sec


    def del_map_spacing_col(self):
        del self.__map_spacing_col


    def del_map_spacing_row(self):
        del self.__map_spacing_row


    def del_map_spacing_sec(self):
        del self.__map_spacing_sec


    def del_map_cell_a(self):
        del self.__map_cell_a


    def del_map_cell_b(self):
        del self.__map_cell_b


    def del_map_cell_c(self):
        del self.__map_cell_c


    def del_map_cell_alpha(self):
        del self.__map_cell_alpha


    def del_map_cell_beta(self):
        del self.__map_cell_beta


    def del_map_cell_gamma(self):
        del self.__map_cell_gamma


    def del_map_axis_order_fast(self):
        del self.__map_axis_order_fast


    def del_map_axis_order_medium(self):
        del self.__map_axis_order_medium


    def del_map_axis_order_slow(self):
        del self.__map_axis_order_slow


    def del_map_minimum(self):
        del self.__map_minimum


    def del_map_maximum(self):
        del self.__map_maximum


    def del_map_average(self):
        del self.__map_average


    def del_map_std(self):
        del self.__map_std


    def del_map_space_group_number(self):
        del self.__map_space_group_number


    def del_map_details(self):
        del self.__map_details


    def del_map_pixel_x(self):
        del self.__map_pixel_x


    def del_map_pixel_y(self):
        del self.__map_pixel_y


    def del_map_pixel_z(self):
        del self.__map_pixel_z


    def del_map_countour_level(self):
        del self.__map_countour_level


    def del_map_annotation_details(self):
        del self.__map_annotation_details

    map_id = property(get_map_id, set_map_id, del_map_id, "map_id's docstring")
    map_file_information = property(get_map_file_information, set_map_file_information, del_map_file_information, "map_file_information's docstring")
    map_data_type = property(get_map_data_type, set_map_data_type, del_map_data_type, "map_data_type's docstring")
    map_num_columns = property(get_map_num_columns, set_map_num_columns, del_map_num_columns, "map_num_columns's docstring")
    map_num_rows = property(get_map_num_rows, set_map_num_rows, del_map_num_rows, "map_num_rows's docstring")
    map_num_sections = property(get_map_num_sections, set_map_num_sections, del_map_num_sections, "map_num_sections's docstring")
    map_origin_col = property(get_map_origin_col, set_map_origin_col, del_map_origin_col, "map_origin_col's docstring")
    map_origin_row = property(get_map_origin_row, set_map_origin_row, del_map_origin_row, "map_origin_row's docstring")
    map_origin_sec = property(get_map_origin_sec, set_map_origin_sec, del_map_origin_sec, "map_origin_sec's docstring")
    map_limit_col = property(get_map_limit_col, set_map_limit_col, del_map_limit_col, "map_limit_col's docstring")
    map_limit_row = property(get_map_limit_row, set_map_limit_row, del_map_limit_row, "map_limit_row's docstring")
    map_limit_sec = property(get_map_limit_sec, set_map_limit_sec, del_map_limit_sec, "map_limit_sec's docstring")
    map_spacing_col = property(get_map_spacing_col, set_map_spacing_col, del_map_spacing_col, "map_spacing_col's docstring")
    map_spacing_row = property(get_map_spacing_row, set_map_spacing_row, del_map_spacing_row, "map_spacing_row's docstring")
    map_spacing_sec = property(get_map_spacing_sec, set_map_spacing_sec, del_map_spacing_sec, "map_spacing_sec's docstring")
    map_cell_a = property(get_map_cell_a, set_map_cell_a, del_map_cell_a, "map_cell_a's docstring")
    map_cell_b = property(get_map_cell_b, set_map_cell_b, del_map_cell_b, "map_cell_b's docstring")
    map_cell_c = property(get_map_cell_c, set_map_cell_c, del_map_cell_c, "map_cell_c's docstring")
    map_cell_alpha = property(get_map_cell_alpha, set_map_cell_alpha, del_map_cell_alpha, "map_cell_alpha's docstring")
    map_cell_beta = property(get_map_cell_beta, set_map_cell_beta, del_map_cell_beta, "map_cell_beta's docstring")
    map_cell_gamma = property(get_map_cell_gamma, set_map_cell_gamma, del_map_cell_gamma, "map_cell_gamma's docstring")
    map_axis_order_fast = property(get_map_axis_order_fast, set_map_axis_order_fast, del_map_axis_order_fast, "map_axis_order_fast's docstring")
    map_axis_order_medium = property(get_map_axis_order_medium, set_map_axis_order_medium, del_map_axis_order_medium, "map_axis_order_medium's docstring")
    map_axis_order_slow = property(get_map_axis_order_slow, set_map_axis_order_slow, del_map_axis_order_slow, "map_axis_order_slow's docstring")
    map_minimum = property(get_map_minimum, set_map_minimum, del_map_minimum, "map_minimum's docstring")
    map_maximum = property(get_map_maximum, set_map_maximum, del_map_maximum, "map_maximum's docstring")
    map_average = property(get_map_average, set_map_average, del_map_average, "map_average's docstring")
    map_std = property(get_map_std, set_map_std, del_map_std, "map_std's docstring")
    map_space_group_number = property(get_map_space_group_number, set_map_space_group_number, del_map_space_group_number, "map_space_group_number's docstring")
    map_details = property(get_map_details, set_map_details, del_map_details, "map_details's docstring")
    map_pixel_x = property(get_map_pixel_x, set_map_pixel_x, del_map_pixel_x, "map_pixel_x's docstring")
    map_pixel_y = property(get_map_pixel_y, set_map_pixel_y, del_map_pixel_y, "map_pixel_y's docstring")
    map_pixel_z = property(get_map_pixel_z, set_map_pixel_z, del_map_pixel_z, "map_pixel_z's docstring")
    map_countour_level = property(get_map_countour_level, set_map_countour_level, del_map_countour_level, "map_countour_level's docstring")
    map_annotation_details = property(get_map_annotation_details, set_map_annotation_details, del_map_annotation_details, "map_annotation_details's docstring")

    
class Emd_entry(Map_information):
    '''
    classdocs
    '''
    __id = None
    __full_name = None
    __acronym = None
    __volume = None
    __resolution = None
    __image_url = None
    __xml_url = None
    __map_url = None

    def __init__(self, id, full_name, acronym, volume, resolution, image_url, xml_url, map_url, map_id, map_file_information, map_data_type, map_num_columns, map_num_rows, map_num_sections, map_origin_col, map_origin_row, map_origin_sec, map_limit_col, map_limit_row, map_limit_sec, map_spacing_col, map_spacing_row, map_spacing_sec, map_cell_a, map_cell_b, map_cell_c, map_cell_alpha, map_cell_beta, map_cell_gamma, map_axis_order_fast, map_axis_order_medium, map_axis_order_slow, map_minimum, map_maximum, map_average, map_std, map_space_group_number, map_details, map_pixel_x, map_pixel_y, map_pixel_z, map_countour_level, map_annotation_details):
        Map_information.__init__(self, 
            map_id, map_file_information, map_data_type, map_num_columns, map_num_rows, map_num_sections, map_origin_col, map_origin_row, map_origin_sec, map_limit_col, map_limit_row, map_limit_sec, map_spacing_col, map_spacing_row, map_spacing_sec, map_cell_a, map_cell_b, map_cell_c, map_cell_alpha, map_cell_beta, map_cell_gamma, map_axis_order_fast, map_axis_order_medium, map_axis_order_slow, map_minimum, map_maximum, map_average, map_std, map_space_group_number, map_details, map_pixel_x, map_pixel_y, map_pixel_z, map_countour_level, map_annotation_details)
        self.__id = id
        self.__full_name = full_name
        self.__acronym = acronym
        self.__volume = volume
        self.__resolution = resolution
        self.__image_url = image_url
        self.__xml_url = xml_url
        self.__map_url = map_url
        
        
    def __insert_map_db(self, cur):
        cur.execute(sql.SQL("INSERT INTO map_information (id,file_information,data_type,num_columns,num_rows,num_sections,origin_col,origin_row,origin_sec,limit_col,limit_row,limit_sec,spacing_col,spacing_row,spacing_sec,cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma,axis_order_fast,axis_order_medium,axis_order_slow,minimum,maximum,average,std,space_group_number,details,pixel_x,pixel_y,pixel_z,countour_level,annotation_details) VALUES (DEFAULT,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s) RETURNING id;")
        ,[json.dumps(self.map_file_information),self.map_data_type,
          self.map_num_columns,self.map_num_rows,self.map_num_sections,
          self.map_origin_col,self.map_origin_row,self.map_origin_sec,
          self.map_limit_col,self.map_limit_row,self.map_limit_sec,
          self.map_spacing_col,self.map_spacing_row,self.map_spacing_sec,
          json.dumps(self.map_cell_a),json.dumps(self.map_cell_b),json.dumps(self.map_cell_c),
          json.dumps(self.map_cell_alpha),json.dumps(self.map_cell_beta),json.dumps(self.map_cell_gamma),
          self.map_axis_order_fast,self.map_axis_order_medium,self.map_axis_order_slow,
          self.map_minimum,self.map_maximum,self.map_average,self.map_std,
          self.map_space_group_number,self.map_details,json.dumps(self.map_pixel_x),
          json.dumps(self.map_pixel_y),json.dumps(self.map_pixel_z),self.map_countour_level,
          self.map_annotation_details])
        self.map_id = [record for record in cur][0]
        
    def __update_map_db(self, cur):
        cur.execute(sql.SQL("UPDATE map_information SET file_information = %s, data_type = %s, num_columns = %s, num_rows = %s, num_sections = %s, origin_col = %s, origin_row = %s, origin_sec = %s, limit_col = %s, limit_row = %s, limit_sec = %s, spacing_col = %s, spacing_row = %s, spacing_sec = %s, cell_a = %s, cell_b = %s, cell_c = %s, cell_alpha = %s, cell_beta = %s, cell_gamma = %s, axis_order_fast = %s, axis_order_medium = %s, axis_order_slow = %s, minimum = %s, maximum = %s, average = %s, std = %s, space_group_number = %s, details = %s, pixel_x = %s, pixel_y = %s, pixel_z = %s, countour_level = %s, annotation_details = %s WHERE id = %s;")
        ,[json.dumps(self.map_file_information),self.map_data_type,
          self.map_num_columns,self.map_num_rows,self.map_num_sections,
          self.map_origin_col,self.map_origin_row,self.map_origin_sec,
          self.map_limit_col,self.map_limit_row,self.map_limit_sec,
          self.map_spacing_col,self.map_spacing_row,self.map_spacing_sec,
          json.dumps(self.map_cell_a),json.dumps(self.map_cell_b),json.dumps(self.map_cell_c),
          json.dumps(self.map_cell_alpha),json.dumps(self.map_cell_beta),json.dumps(self.map_cell_gamma),
          self.map_axis_order_fast,self.map_axis_order_medium,self.map_axis_order_slow,
          self.map_minimum,self.map_maximum,self.map_average,self.map_std,
          self.map_space_group_number,self.map_details,json.dumps(self.map_pixel_x),
          json.dumps(self.map_pixel_y),json.dumps(self.map_pixel_z),self.map_countour_level,
          self.map_annotation_details,self.map_id])
    
    def __insert_emd_entry(self, cur):
        cur.execute(sql.SQL("INSERT INTO emd_entry(id,full_name,acronym,volume,resolution,image_url,xml_url,map_url,map_information_id) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s);")
        ,[self.id,self.full_name,self.acronym,self.volume,self.resolution,self.image_url,self.xml_url,self.map_url,self.map_id])
        
    def __update_emd_entry(self, cur):
        cur.execute(sql.SQL("UPDATE emd_entry SET full_name = %s, acronym = %s, volume = %s, resolution = %s, image_url = %s, xml_url = %s, map_url = %s, map_information_id = %s WHERE id = %s;")
        ,[self.full_name,self.acronym,self.volume,self.resolution,self.image_url,self.xml_url,self.map_url,self.map_id,self.id])    
        
    def insert_db(self, cur):
        self.__insert_map_db(cur)
        self.__insert_emd_entry(cur)
        
    def update_db(self, cur):
        self.__update_map_db(cur)
        self.__update_emd_entry(cur)

    def get_id(self):
        return self.__id


    def get_full_name(self):
        return self.__full_name


    def get_acronym(self):
        return self.__acronym


    def get_volume(self):
        return self.__volume


    def get_resolution(self):
        return self.__resolution


    def get_image_url(self):
        return self.__image_url


    def get_xml_url(self):
        return self.__xml_url


    def get_map_url(self):
        return self.__map_url


    def set_id(self, value):
        self.__id = value


    def set_full_name(self, value):
        self.__full_name = value


    def set_acronym(self, value):
        self.__acronym = value


    def set_volume(self, value):
        self.__volume = value


    def set_resolution(self, value):
        self.__resolution = value


    def set_image_url(self, value):
        self.__image_url = value


    def set_xml_url(self, value):
        self.__xml_url = value


    def set_map_url(self, value):
        self.__map_url = value


    def del_id(self):
        del self.__id


    def del_full_name(self):
        del self.__full_name


    def del_acronym(self):
        del self.__acronym


    def del_volume(self):
        del self.__volume


    def del_resolution(self):
        del self.__resolution


    def del_image_url(self):
        del self.__image_url


    def del_xml_url(self):
        del self.__xml_url


    def del_map_url(self):
        del self.__map_url
        
        
    def __eq__(self, emd_entry):        
        return self.id == emd_entry.id \
           and self.full_name == emd_entry.full_name \
           and self.acronym == emd_entry.acronym \
           and self.volume == emd_entry.volume \
           and self.resolution == emd_entry.resolution \
           and self.image_url == emd_entry.image_url \
           and self.xml_url == emd_entry.xml_url \
           and self.map_url == emd_entry.map_url \
           and self.map_id == emd_entry.map_id \
           and self.map_file_information == emd_entry.map_file_information \
           and self.map_data_type == emd_entry.map_data_type \
           and self.map_num_columns == emd_entry.map_num_columns \
           and self.map_num_rows == emd_entry.map_num_rows \
           and self.map_num_sections == emd_entry.map_num_sections \
           and self.map_origin_col == emd_entry.map_origin_col \
           and self.map_origin_row == emd_entry.map_origin_row \
           and self.map_origin_sec == emd_entry.map_origin_sec \
           and self.map_limit_col == emd_entry.map_limit_col \
           and self.map_limit_row == emd_entry.map_limit_row \
           and self.map_limit_sec == emd_entry.map_limit_sec \
           and self.map_spacing_col == emd_entry.map_spacing_col \
           and self.map_spacing_row == emd_entry.map_spacing_row \
           and self.map_spacing_sec == emd_entry.map_spacing_sec \
           and self.map_cell_a == emd_entry.map_cell_a \
           and self.map_cell_b == emd_entry.map_cell_b \
           and self.map_cell_c == emd_entry.map_cell_c \
           and self.map_cell_alpha == emd_entry.map_cell_alpha \
           and self.map_cell_beta == emd_entry.map_cell_beta \
           and self.map_cell_gamma == emd_entry.map_cell_gamma \
           and self.map_axis_order_fast == emd_entry.map_axis_order_fast \
           and self.map_axis_order_medium == emd_entry.map_axis_order_medium \
           and self.map_axis_order_slow == emd_entry.map_axis_order_slow \
           and self.map_minimum == emd_entry.map_minimum \
           and self.map_maximum == emd_entry.map_maximum \
           and self.map_average == emd_entry.map_average \
           and self.map_std == emd_entry.map_std \
           and self.map_space_group_number == emd_entry.map_space_group_number \
           and self.map_details == emd_entry.map_details \
           and self.map_pixel_x == emd_entry.map_pixel_x \
           and self.map_pixel_y == emd_entry.map_pixel_y \
           and self.map_pixel_z == emd_entry.map_pixel_z \
           and self.map_countour_level == emd_entry.map_countour_level \
           and self.map_annotation_details == emd_entry.map_annotation_details

    id = property(get_id, set_id, del_id, "id's docstring")
    full_name = property(get_full_name, set_full_name, del_full_name, "full_name's docstring")
    acronym = property(get_acronym, set_acronym, del_acronym, "acronym's docstring")
    volume = property(get_volume, set_volume, del_volume, "volume's docstring")
    resolution = property(get_resolution, set_resolution, del_resolution, "resolution's docstring")
    image_url = property(get_image_url, set_image_url, del_image_url, "image_url's docstring")
    xml_url = property(get_xml_url, set_xml_url, del_xml_url, "xml_url's docstring")
    map_url = property(get_map_url, set_map_url, del_map_url, "map_url's docstring")

    
    