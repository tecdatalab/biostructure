'''
Created on 5 oct. 2019

@author: luis98
'''
from psycopg2 import sql


class Cath_atomic_structure(object):
    '''
    classdocs
    ''' 
    __atomic_structure_id = None
    __id_code = None
    __class_number = None
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
    
    def __init__(self, atomic_structure_id = None, id_code = None, class_number = None, architecture_number = None, topology_number = None, homologous_superfamily_number = None, S35_sequence_cluster_number = None, S60_sequence_cluster_number = None, S95_sequence_cluster_number = None, S100_sequence_cluster_number = None, S100_sequence_count_number = None, domain_length = None, structure_resolution = None):
        self.__atomic_structure_id = atomic_structure_id
        self.__id_code = id_code
        self.__class_number = class_number
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
        
    def insert_db(self, cur):
        if self.__atomic_structure_id == None:
            cur.execute(
            sql.SQL("SELECT id FROM atomic_structure WHERE id_code = %s;"),[
                self.__id_code])
            result = [record[0] for record in cur]
            self.__atomic_structure_id= result[0]
        
        cur.execute(
            sql.SQL("INSERT INTO cath_atomic_structure(atomic_structure_id,class_number,architecture_number,topology_number,homologous_superfamily_number,S35_sequence_cluster_number,S60_sequence_cluster_number,S95_sequence_cluster_number,S100_sequence_cluster_number,S100_sequence_count_number,domain_length,structure_resolution) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s); "),[
                self.__atomic_structure_id,
                self.__class_number,
                self.__architecture_number,
                self.__topology_number,
                self.__homologous_superfamily_number,
                self.__S35_sequence_cluster_number,
                self.__S60_sequence_cluster_number,
                self.__S95_sequence_cluster_number,
                self.__S100_sequence_cluster_number,
                self.__S100_sequence_count_number,
                self.__domain_length,
                self.__structure_resolution])
        
    def update_db(self, cur):
        if self.atomic_structure_id == None:
            cur.execute(
            sql.SQL("SELECT id FROM atomic_structure WHERE id_code = %s;"),[
                self.__id_code])
            result = [record[0] for record in cur]
            self.__atomic_structure_id= result[0]
        
        cur.execute(
            sql.SQL("UPDATE cath_atomic_structure SET class_number = %s,architecture_number = %s,topology_number = %s,homologous_superfamily_number = %s,S35_sequence_cluster_number = %s,S60_sequence_cluster_number = %s,S95_sequence_cluster_number = %s,S100_sequence_cluster_number = %s,S100_sequence_count_number = %s,domain_length = %s,structure_resolution = %s WHERE atomic_structure_id = %s; "),[
                self.__class_number,
                self.__architecture_number,
                self.__topology_number,
                self.__homologous_superfamily_number,
                self.__S35_sequence_cluster_number,
                self.__S60_sequence_cluster_number,
                self.__S95_sequence_cluster_number,
                self.__S100_sequence_cluster_number,
                self.__S100_sequence_count_number,
                self.__domain_length,
                self.__structure_resolution,
                self.__atomic_structure_id])
    
    def insert_update_db(self, cur):
        if self.__atomic_structure_id == None:
            cur.execute(
            sql.SQL("SELECT id FROM atomic_structure WHERE id_code = %s;"),[
                self.__id_code])
            result = [record[0] for record in cur]
            self.__atomic_structure_id= result[0]
        
        cur.execute(
            sql.SQL("SELECT * FROM cath_atomic_structure WHERE atomic_structure_id = %s;"),[
                self.__atomic_structure_id])
        result = [record[0] for record in cur]
 
        if len(result)>0:
            self.update_db(cur)
        else:
            self.insert_db(cur)
        
    def get_id_code(self):
        return self.__id_code


    def get_class_number(self):
        return self.__class_number


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


    def get_atomic_structure_id(self):
        return self.__atomic_structure_id


    def set_id_code(self, value):
        self.__id_code = value


    def set_class_number(self, value):
        self.__class_number = value


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


    def set_atomic_structure_id(self, value):
        self.__atomic_structure_id = value



    def del_id_code(self):
        del self.__id_code


    def del_class_number(self):
        del self.__class_number


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


    def del_atomic_structure_id(self):
        del self.__atomic_structure_id

    id_code = property(get_id_code, set_id_code, del_id_code, "id_code's docstring")
    class_number = property(get_class_number, set_class_number, del_class_number, "class_number's docstring")
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
    atomic_structure_id = property(get_atomic_structure_id, set_atomic_structure_id, del_atomic_structure_id, "atomic_structure_id's docstring")




    