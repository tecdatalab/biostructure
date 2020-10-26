import sys
from em.dataset.feature_extraction import annotateSample
from em.molecule import Molecule 
import pandas as pd


def test_annotate_segments():
    # Get id unique values to extract indexes of respective molecule subunits 
    entries = [['maps/segments.map','1','maps/segments_A.map','1'],
               ['maps/segments.map','1','maps/segments_B.map','2'],
               ['maps/segments.map','1','maps/segments_C.map','3'],
               ['maps/segments.map','1','maps/segments_D.map','4'],
               ['maps/segments.map','1','maps/segments_E.map','5'],
               ['maps/segments.map','1','maps/segments_F.map','6'],
               ['maps/segments.map','1','maps/segments_G.map','7'],
               ['maps/segments.map','1','maps/segments_H.map','8']]
    df = pd.DataFrame(entries,columns=['map_path','contourLevel','subunit_path','chain_label'])
    indexes = df.index.tolist()
    result = annotateSample(indexes,df,3,'./')
    
    structure = Molecule('maps/segments.map',1)
    structure_vol = structure.getVolume()[1]
 
    segment = Molecule('maps/segments_A.map',1)
    expected_volume = segment.getVolume()[1]
    assert result=={'total':structure_vol, 1:expected_volume, 2:expected_volume, 3:expected_volume, 4:expected_volume,
                                           5:expected_volume, 6:expected_volume, 7:expected_volume, 8:expected_volume}
    
    
