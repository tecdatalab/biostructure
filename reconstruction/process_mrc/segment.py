'''
Created on Jun 9, 2020

@author: luis98
'''


class Segment(object):
  '''
  classdocs
  '''

  def __init__(self, id_segment, mask, zd_descriptors, volume=0, textSimulatePDB=None, recommendedContour_p=0):
    '''
    Constructor
    '''
    self.id_segment = id_segment
    self.mask = mask
    self.zd_descriptors = zd_descriptors
    self.volume = volume
    self.textSimulatePDB = textSimulatePDB
    self.recommendedContour = recommendedContour_p
    self.originalCA = None
    self.simulateCA = None
