'''
Created on 22 feb. 2019

@author: luis98
'''

class Map_information(object):
    '''
    classdocs
    '''
    __id = None
    __file_information = None
    __data_type = None
    __num_columns = None
    __num_rows = None
    __num_sections = None
    __origin_col = None
    __origin_row = None
    __origin_sec = None
    __limit_col = None
    __limit_row = None
    __limit_sec = None
    __spacing_col = None
    __spacing_row = None
    __spacing_sec = None
    __cell_a = None
    __cell_b = None
    __cell_c = None
    __cell_alpha = None
    __cell_beta = None
    __cell_gamma = None
    __axis_order_fast = None
    __axis_order_medium = None
    __axis_order_slow = None
    __minimum = None
    __maximum = None
    __average = None
    __std = None
    __space_group_number = None
    __datails = None
    __pixel_x = None
    __pixel_y = None
    __pixel_z = None
    __countour_level = None
    __annotation_details = None

    def __init__(self, id, file_information, data_type, num_columns, num_rows, num_sections, origin_col, origin_row, origin_sec, limit_col, limit_row, limit_sec, spacing_col, spacing_row, spacing_sec, cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma, axis_order_fast, axis_order_medium, axis_order_slow, minimum, maximum, average, std, space_group_number, datails, pixel_x, pixel_y, pixel_z, countour_level, annotation_details):
        self.__id = id
        self.__file_information = file_information
        self.__data_type = data_type
        self.__num_columns = num_columns
        self.__num_rows = num_rows
        self.__num_sections = num_sections
        self.__origin_col = origin_col
        self.__origin_row = origin_row
        self.__origin_sec = origin_sec
        self.__limit_col = limit_col
        self.__limit_row = limit_row
        self.__limit_sec = limit_sec
        self.__spacing_col = spacing_col
        self.__spacing_row = spacing_row
        self.__spacing_sec = spacing_sec
        self.__cell_a = cell_a
        self.__cell_b = cell_b
        self.__cell_c = cell_c
        self.__cell_alpha = cell_alpha
        self.__cell_beta = cell_beta
        self.__cell_gamma = cell_gamma
        self.__axis_order_fast = axis_order_fast
        self.__axis_order_medium = axis_order_medium
        self.__axis_order_slow = axis_order_slow
        self.__minimum = minimum
        self.__maximum = maximum
        self.__average = average
        self.__std = std
        self.__space_group_number = space_group_number
        self.__datails = datails
        self.__pixel_x = pixel_x
        self.__pixel_y = pixel_y
        self.__pixel_z = pixel_z
        self.__countour_level = countour_level
        self.__annotation_details = annotation_details

    def get_id(self):
        return self.__id


    def get_file_information(self):
        return self.__file_information


    def get_data_type(self):
        return self.__data_type


    def get_num_columns(self):
        return self.__num_columns


    def get_num_rows(self):
        return self.__num_rows


    def get_num_sections(self):
        return self.__num_sections


    def get_origin_col(self):
        return self.__origin_col


    def get_origin_row(self):
        return self.__origin_row


    def get_origin_sec(self):
        return self.__origin_sec


    def get_limit_col(self):
        return self.__limit_col


    def get_limit_row(self):
        return self.__limit_row


    def get_limit_sec(self):
        return self.__limit_sec


    def get_spacing_col(self):
        return self.__spacing_col


    def get_spacing_row(self):
        return self.__spacing_row


    def get_spacing_sec(self):
        return self.__spacing_sec


    def get_cell_a(self):
        return self.__cell_a


    def get_cell_b(self):
        return self.__cell_b


    def get_cell_c(self):
        return self.__cell_c


    def get_cell_alpha(self):
        return self.__cell_alpha


    def get_cell_beta(self):
        return self.__cell_beta


    def get_cell_gamma(self):
        return self.__cell_gamma


    def get_axis_order_fast(self):
        return self.__axis_order_fast


    def get_axis_order_medium(self):
        return self.__axis_order_medium


    def get_axis_order_slow(self):
        return self.__axis_order_slow


    def get_minimum(self):
        return self.__minimum


    def get_maximum(self):
        return self.__maximum


    def get_average(self):
        return self.__average


    def get_std(self):
        return self.__std


    def get_space_group_number(self):
        return self.__space_group_number


    def get_datails(self):
        return self.__datails


    def get_pixel_x(self):
        return self.__pixel_x


    def get_pixel_y(self):
        return self.__pixel_y


    def get_pixel_z(self):
        return self.__pixel_z


    def get_countour_level(self):
        return self.__countour_level


    def get_annotation_details(self):
        return self.__annotation_details


    def set_id(self, value):
        self.__id = value


    def set_file_information(self, value):
        self.__file_information = value


    def set_data_type(self, value):
        self.__data_type = value


    def set_num_columns(self, value):
        self.__num_columns = value


    def set_num_rows(self, value):
        self.__num_rows = value


    def set_num_sections(self, value):
        self.__num_sections = value


    def set_origin_col(self, value):
        self.__origin_col = value


    def set_origin_row(self, value):
        self.__origin_row = value


    def set_origin_sec(self, value):
        self.__origin_sec = value


    def set_limit_col(self, value):
        self.__limit_col = value


    def set_limit_row(self, value):
        self.__limit_row = value


    def set_limit_sec(self, value):
        self.__limit_sec = value


    def set_spacing_col(self, value):
        self.__spacing_col = value


    def set_spacing_row(self, value):
        self.__spacing_row = value


    def set_spacing_sec(self, value):
        self.__spacing_sec = value


    def set_cell_a(self, value):
        self.__cell_a = value


    def set_cell_b(self, value):
        self.__cell_b = value


    def set_cell_c(self, value):
        self.__cell_c = value


    def set_cell_alpha(self, value):
        self.__cell_alpha = value


    def set_cell_beta(self, value):
        self.__cell_beta = value


    def set_cell_gamma(self, value):
        self.__cell_gamma = value


    def set_axis_order_fast(self, value):
        self.__axis_order_fast = value


    def set_axis_order_medium(self, value):
        self.__axis_order_medium = value


    def set_axis_order_slow(self, value):
        self.__axis_order_slow = value


    def set_minimum(self, value):
        self.__minimum = value


    def set_maximum(self, value):
        self.__maximum = value


    def set_average(self, value):
        self.__average = value


    def set_std(self, value):
        self.__std = value


    def set_space_group_number(self, value):
        self.__space_group_number = value


    def set_datails(self, value):
        self.__datails = value


    def set_pixel_x(self, value):
        self.__pixel_x = value


    def set_pixel_y(self, value):
        self.__pixel_y = value


    def set_pixel_z(self, value):
        self.__pixel_z = value


    def set_countour_level(self, value):
        self.__countour_level = value


    def set_annotation_details(self, value):
        self.__annotation_details = value


    def del_id(self):
        del self.__id


    def del_file_information(self):
        del self.__file_information


    def del_data_type(self):
        del self.__data_type


    def del_num_columns(self):
        del self.__num_columns


    def del_num_rows(self):
        del self.__num_rows


    def del_num_sections(self):
        del self.__num_sections


    def del_origin_col(self):
        del self.__origin_col


    def del_origin_row(self):
        del self.__origin_row


    def del_origin_sec(self):
        del self.__origin_sec


    def del_limit_col(self):
        del self.__limit_col


    def del_limit_row(self):
        del self.__limit_row


    def del_limit_sec(self):
        del self.__limit_sec


    def del_spacing_col(self):
        del self.__spacing_col


    def del_spacing_row(self):
        del self.__spacing_row


    def del_spacing_sec(self):
        del self.__spacing_sec


    def del_cell_a(self):
        del self.__cell_a


    def del_cell_b(self):
        del self.__cell_b


    def del_cell_c(self):
        del self.__cell_c


    def del_cell_alpha(self):
        del self.__cell_alpha


    def del_cell_beta(self):
        del self.__cell_beta


    def del_cell_gamma(self):
        del self.__cell_gamma


    def del_axis_order_fast(self):
        del self.__axis_order_fast


    def del_axis_order_medium(self):
        del self.__axis_order_medium


    def del_axis_order_slow(self):
        del self.__axis_order_slow


    def del_minimum(self):
        del self.__minimum


    def del_maximum(self):
        del self.__maximum


    def del_average(self):
        del self.__average


    def del_std(self):
        del self.__std


    def del_space_group_number(self):
        del self.__space_group_number


    def del_datails(self):
        del self.__datails


    def del_pixel_x(self):
        del self.__pixel_x


    def del_pixel_y(self):
        del self.__pixel_y


    def del_pixel_z(self):
        del self.__pixel_z


    def del_countour_level(self):
        del self.__countour_level


    def del_annotation_details(self):
        del self.__annotation_details

    id = property(get_id, set_id, del_id, "id's docstring")
    file_information = property(get_file_information, set_file_information, del_file_information, "file_information's docstring")
    data_type = property(get_data_type, set_data_type, del_data_type, "data_type's docstring")
    num_columns = property(get_num_columns, set_num_columns, del_num_columns, "num_columns's docstring")
    num_rows = property(get_num_rows, set_num_rows, del_num_rows, "num_rows's docstring")
    num_sections = property(get_num_sections, set_num_sections, del_num_sections, "num_sections's docstring")
    origin_col = property(get_origin_col, set_origin_col, del_origin_col, "origin_col's docstring")
    origin_row = property(get_origin_row, set_origin_row, del_origin_row, "origin_row's docstring")
    origin_sec = property(get_origin_sec, set_origin_sec, del_origin_sec, "origin_sec's docstring")
    limit_col = property(get_limit_col, set_limit_col, del_limit_col, "limit_col's docstring")
    limit_row = property(get_limit_row, set_limit_row, del_limit_row, "limit_row's docstring")
    limit_sec = property(get_limit_sec, set_limit_sec, del_limit_sec, "limit_sec's docstring")
    spacing_col = property(get_spacing_col, set_spacing_col, del_spacing_col, "spacing_col's docstring")
    spacing_row = property(get_spacing_row, set_spacing_row, del_spacing_row, "spacing_row's docstring")
    spacing_sec = property(get_spacing_sec, set_spacing_sec, del_spacing_sec, "spacing_sec's docstring")
    cell_a = property(get_cell_a, set_cell_a, del_cell_a, "cell_a's docstring")
    cell_b = property(get_cell_b, set_cell_b, del_cell_b, "cell_b's docstring")
    cell_c = property(get_cell_c, set_cell_c, del_cell_c, "cell_c's docstring")
    cell_alpha = property(get_cell_alpha, set_cell_alpha, del_cell_alpha, "cell_alpha's docstring")
    cell_beta = property(get_cell_beta, set_cell_beta, del_cell_beta, "cell_beta's docstring")
    cell_gamma = property(get_cell_gamma, set_cell_gamma, del_cell_gamma, "cell_gamma's docstring")
    axis_order_fast = property(get_axis_order_fast, set_axis_order_fast, del_axis_order_fast, "axis_order_fast's docstring")
    axis_order_medium = property(get_axis_order_medium, set_axis_order_medium, del_axis_order_medium, "axis_order_medium's docstring")
    axis_order_slow = property(get_axis_order_slow, set_axis_order_slow, del_axis_order_slow, "axis_order_slow's docstring")
    minimum = property(get_minimum, set_minimum, del_minimum, "minimum's docstring")
    maximum = property(get_maximum, set_maximum, del_maximum, "maximum's docstring")
    average = property(get_average, set_average, del_average, "average's docstring")
    std = property(get_std, set_std, del_std, "std's docstring")
    space_group_number = property(get_space_group_number, set_space_group_number, del_space_group_number, "space_group_number's docstring")
    datails = property(get_datails, set_datails, del_datails, "datails's docstring")
    pixel_x = property(get_pixel_x, set_pixel_x, del_pixel_x, "pixel_x's docstring")
    pixel_y = property(get_pixel_y, set_pixel_y, del_pixel_y, "pixel_y's docstring")
    pixel_z = property(get_pixel_z, set_pixel_z, del_pixel_z, "pixel_z's docstring")
    countour_level = property(get_countour_level, set_countour_level, del_countour_level, "countour_level's docstring")
    annotation_details = property(get_annotation_details, set_annotation_details, del_annotation_details, "annotation_details's docstring")

    
    