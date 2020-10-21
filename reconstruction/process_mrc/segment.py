'''
Created on Jun 9, 2020

@author: luis98
'''


class Segment(object):
    '''
    classdocs
    '''

    def __init__(self, id_segment, mask, zd_descriptors):
        '''
        Constructor
        '''
        self.id_segment = id_segment
        self.mask = mask
        self.zd_descriptors = zd_descriptors
