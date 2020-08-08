'''
Created on 7 aug. 2020

@author: dnnxl
'''
from psycopg2 import sql

class Algorithm_entry(object):
    '''
    The Algorithm_entry object deal with the algorithm to generate the segments.

    Parameters:
        None
    '''
    def get_algorithm_id(self, cur, algorithm):
        '''
        Returns the id to the correspond algorithm name.

        Parameters:
            cur (object)         -  The sql cursor object.
            algorithm (string)   -  The algorithm name.
        Returns:
            result[0](int or none)  -  The id of the algorithm name.
        '''
        cur.execute(
            sql.SQL("SELECT id_algorithm FROM segment_algorithm WHERE algorithm_name = %s;"), [algorithm])
        result = [record[0] for record in cur]
        if len(result)>0:
            return result[0]
        else:
            return None