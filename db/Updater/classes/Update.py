'''
Created on 22 feb. 2019

@author: luis98
'''

class Update(object):
    '''
    classdocs
    '''
    __last_update = None

    def __init__(self, last_update):
        self.__last_update = last_update

    def get_last_update(self):
        return self.__last_update


    def set_last_update(self, value):
        self.__last_update = value


    def del_last_update(self):
        del self.__last_update

    last_update = property(get_last_update, set_last_update, del_last_update, "last_update's docstring")      