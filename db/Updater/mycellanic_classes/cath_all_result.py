'''
Created on 5 oct. 2019

@author: luis98
'''

class Cath_all_result(object):
    '''
    classdocs
    ''' 
    __CATH_domain_name = None
    __lass_number = None
    __architecture_number = None
    __topology_number = None
    __homologous_superfamily_number = None
    __S35_sequence_cluster_number = None
    __S60_sequence_cluster_number = None
    __S95_sequence_cluster_number = None
    __S100_sequence_cluster_number = None
    __S100_sequence_count_number = None
    __domain_length = None
    __structure_resolution = None
    
    def __init__(self, CATH_domain_name, lass_number, architecture_number, topology_number, homologous_superfamily_number, S35_sequence_cluster_number, S60_sequence_cluster_number, S95_sequence_cluster_number, S100_sequence_cluster_number, S100_sequence_count_number, domain_length, structure_resolution):
        self.__CATH_domain_name = CATH_domain_name
        self.__lass_number = lass_number
        self.__architecture_number = architecture_number
        self.__topology_number = topology_number
        self.__homologous_superfamily_number = homologous_superfamily_number
        self.__S35_sequence_cluster_number = S35_sequence_cluster_number
        self.__S60_sequence_cluster_number = S60_sequence_cluster_number
        self.__S95_sequence_cluster_number = S95_sequence_cluster_number
        self.__S100_sequence_cluster_number = S100_sequence_cluster_number
        self.__S100_sequence_count_number = S100_sequence_count_number
        self.__domain_length = domain_length
        self.__structure_resolution = structure_resolution

    def get_cath_domain_name(self):
        return self.__CATH_domain_name


    def get_lass_number(self):
        return self.__lass_number


    def get_architecture_number(self):
        return self.__architecture_number


    def get_topology_number(self):
        return self.__topology_number


    def get_homologous_superfamily_number(self):
        return self.__homologous_superfamily_number


    def get_s_35_sequence_cluster_number(self):
        return self.__S35_sequence_cluster_number


    def get_s_60_sequence_cluster_number(self):
        return self.__S60_sequence_cluster_number


    def get_s_95_sequence_cluster_number(self):
        return self.__S95_sequence_cluster_number


    def get_s_100_sequence_cluster_number(self):
        return self.__S100_sequence_cluster_number


    def get_s_100_sequence_count_number(self):
        return self.__S100_sequence_count_number


    def get_domain_length(self):
        return self.__domain_length


    def get_structure_resolution(self):
        return self.__structure_resolution


    def set_cath_domain_name(self, value):
        self.__CATH_domain_name = value


    def set_lass_number(self, value):
        self.__lass_number = value


    def set_architecture_number(self, value):
        self.__architecture_number = value


    def set_topology_number(self, value):
        self.__topology_number = value


    def set_homologous_superfamily_number(self, value):
        self.__homologous_superfamily_number = value


    def set_s_35_sequence_cluster_number(self, value):
        self.__S35_sequence_cluster_number = value


    def set_s_60_sequence_cluster_number(self, value):
        self.__S60_sequence_cluster_number = value


    def set_s_95_sequence_cluster_number(self, value):
        self.__S95_sequence_cluster_number = value


    def set_s_100_sequence_cluster_number(self, value):
        self.__S100_sequence_cluster_number = value


    def set_s_100_sequence_count_number(self, value):
        self.__S100_sequence_count_number = value


    def set_domain_length(self, value):
        self.__domain_length = value


    def set_structure_resolution(self, value):
        self.__structure_resolution = value


    def del_cath_domain_name(self):
        del self.__CATH_domain_name


    def del_lass_number(self):
        del self.__lass_number


    def del_architecture_number(self):
        del self.__architecture_number


    def del_topology_number(self):
        del self.__topology_number


    def del_homologous_superfamily_number(self):
        del self.__homologous_superfamily_number


    def del_s_35_sequence_cluster_number(self):
        del self.__S35_sequence_cluster_number


    def del_s_60_sequence_cluster_number(self):
        del self.__S60_sequence_cluster_number


    def del_s_95_sequence_cluster_number(self):
        del self.__S95_sequence_cluster_number


    def del_s_100_sequence_cluster_number(self):
        del self.__S100_sequence_cluster_number


    def del_s_100_sequence_count_number(self):
        del self.__S100_sequence_count_number


    def del_domain_length(self):
        del self.__domain_length


    def del_structure_resolution(self):
        del self.__structure_resolution

    CATH_domain_name = property(get_cath_domain_name, set_cath_domain_name, del_cath_domain_name, "CATH_domain_name's docstring")
    lass_number = property(get_lass_number, set_lass_number, del_lass_number, "lass_number's docstring")
    architecture_number = property(get_architecture_number, set_architecture_number, del_architecture_number, "architecture_number's docstring")
    topology_number = property(get_topology_number, set_topology_number, del_topology_number, "topology_number's docstring")
    homologous_superfamily_number = property(get_homologous_superfamily_number, set_homologous_superfamily_number, del_homologous_superfamily_number, "homologous_superfamily_number's docstring")
    S35_sequence_cluster_number = property(get_s_35_sequence_cluster_number, set_s_35_sequence_cluster_number, del_s_35_sequence_cluster_number, "S35_sequence_cluster_number's docstring")
    S60_sequence_cluster_number = property(get_s_60_sequence_cluster_number, set_s_60_sequence_cluster_number, del_s_60_sequence_cluster_number, "S60_sequence_cluster_number's docstring")
    S95_sequence_cluster_number = property(get_s_95_sequence_cluster_number, set_s_95_sequence_cluster_number, del_s_95_sequence_cluster_number, "S95_sequence_cluster_number's docstring")
    S100_sequence_cluster_number = property(get_s_100_sequence_cluster_number, set_s_100_sequence_cluster_number, del_s_100_sequence_cluster_number, "S100_sequence_cluster_number's docstring")
    S100_sequence_count_number = property(get_s_100_sequence_count_number, set_s_100_sequence_count_number, del_s_100_sequence_count_number, "S100_sequence_count_number's docstring")
    domain_length = property(get_domain_length, set_domain_length, del_domain_length, "domain_length's docstring")
    structure_resolution = property(get_structure_resolution, set_structure_resolution, del_structure_resolution, "structure_resolution's docstring")

