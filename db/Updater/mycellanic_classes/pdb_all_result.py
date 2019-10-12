'''
Created on 22 sep. 2019

@author: luis98
'''


class Pdb_all_result:
    '''
    classdocs
    '''
    __data1 = None
    __data2 = None
    __data3 = None

    def __init__(self, data1, data2, data3):
        self.__data1 = data1
        self.__data2 = data2
        self.__data3 = data3

    def get_data_1(self):
        return self.__data1

    def get_data_2(self):
        return self.__data2

    def get_data_3(self):
        return self.__data3

    def set_data_1(self, value):
        self.__data1 = value

    def set_data_2(self, value):
        self.__data2 = value

    def set_data_3(self, value):
        self.__data3 = value

    def del_data_1(self):
        del self.__data1

    def del_data_2(self):
        del self.__data2

    def del_data_3(self):
        del self.__data3

    data1 = property(get_data_1, set_data_1, del_data_1, "data1's docstring")
    data2 = property(get_data_2, set_data_2, del_data_2, "data2's docstring")
    data3 = property(get_data_3, set_data_3, del_data_3, "data3's docstring")
