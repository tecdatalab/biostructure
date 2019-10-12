'''
Created on 5 oct. 2019

@author: luis98
'''


class Cath_domain_boundaries_result(object):
    __domain = None
    __sequences = []

    def __init__(self, domain, sequences):
        self.__domain = domain
        self.__sequences = sequences

    def get_domain(self):
        return self.__domain

    def get_sequences(self):
        return self.__sequences

    def set_domain(self, value):
        self.__domain = value

    def set_sequences(self, value):
        self.__sequences = value

    def del_domain(self):
        del self.__domain

    def del_sequences(self):
        del self.__sequences

    domain = property(get_domain, set_domain, del_domain, "domain's docstring")
    sequences = property(
        get_sequences,
        set_sequences,
        del_sequences,
        "sequences's docstring")
