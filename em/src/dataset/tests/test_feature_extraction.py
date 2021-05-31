import sys
from em.dataset.feature_extraction import annotateSample
from em.molecule import Molecule 
import pandas as pd
import random

def test_annotate_segments():
    # Get id unique values to extract indexes of respective molecule subunits 
    entries = [['map01','maps/segments.map','1','maps/segments_A.map','1'],
               ['map01','maps/segments.map','1','maps/segments_B.map','2'],
               ['map01','maps/segments.map','1','maps/segments_C.map','3'],
               ['map01','maps/segments.map','1','maps/segments_D.map','4'],
               ['map01','maps/segments.map','1','maps/segments_E.map','5'],
               ['map01','maps/segments.map','1','maps/segments_F.map','6'],
               ['map01','maps/segments.map','1','maps/segments_G.map','7'],
               ['map01','maps/segments.map','1','maps/segments_H.map','8']]
    df = pd.DataFrame(entries,columns=['id','aligned_path','contourLevel','subunit_path','chain_label'])
    indexes = df.index.tolist()
    result = annotateSample('map01', indexes,df,3,{'id':'id','map_path':'aligned_path','contourLevel':'contourLevel', 'subunit_path':'subunit_path','chain_label':'chain_label'}, './')
    
    structure = Molecule('maps/segments.map',1)
    structure_vol = structure.getVolume()[1]
 
    segment = Molecule('maps/segments_A.map',1)
    expected_volume = segment.getVolume()[1]
    json = '1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}, 7: {}, 8: {}'.format(str(expected_volume), str(expected_volume), str(expected_volume), str(expected_volume), str(expected_volume), str(expected_volume), str(expected_volume), str(expected_volume))

    assert result=={'map_id': 'map01','total': 512000, 'voxels_descarted':0, 'voxels_reasigned':0, 'voxels_assigned': '{'+json+'}', 'tag_path':'./maps/segments_gt.map'}
   
'''
def test_annotate_real_data():
    # load csv and sample one map
    df_full = pd.read_csv('/work/mzumbado/biostructure/em/dataset/dataset_selected.csv')
    map_path_list = df_full.aligned_path.unique().tolist()
    map_path = random.choice(map_path_list)
    map_path = map_path_list[0]
    entries = df_full[df_full['aligned_path'] == map_path]
    df = entries[['id','aligned_path','contourLevel','subunit_path','chain_label']]
    indexes = df.index.tolist()
    result = annotateSample(''indexes,df,3,'./',{'id':'id','map_path':'aligned_path','contourLevel':'contourLevel', 'subunit_path':'subunit_path','chain_label':'chain_label'},'./')
    print(result) 
    assert False
'''        
    
