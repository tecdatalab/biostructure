import sys 
sys.path.append('..')
import molecule
import numpy as np

# Segmentation test with 4 regions
def test_segmentation_4_regions():
    segment_map = molecule.Molecule('maps/segments.map', 1)
    segment_map.generateSegments(3,1)
    masks = segment_map.getSegmentsMasks()[1]
    assert len(masks) == 8
    print("test_segmentation_4_regions OK")

# Segmentation test with 1 region
def test_segmentation_1_region():
    map_single_region = molecule.Molecule('maps/cube100x100x100.map', 1)
    map_single_region.generateSegments(3,1)
    masks = map_single_region.getSegmentsMasks()[1]
    assert len(masks) == 1
    print("test_segmentation_1_region OK")

# Test volume calculation
def test_volume_whole():
    map_100 = molecule.Molecule('maps/cube100x100x100.map', 1)
    vol = map_100.getVolume()[1]
    expected = np.sum(map_100.emMap.data() >= 1)
    assert vol == expected
    print('test_volume_whole OK')

# Test volumes at  diferents isosurface level
def test_volume_isosurfaces():

    map_levels = molecule.Molecule('maps/cubeincube.map', 3, [1, 0.5, 0.3])
    vol_levels = map_levels.getVolume()
    expected_at_1 = np.sum((map_levels.emMap.data() == 3)| (map_levels.emMap.data() == 2) | (map_levels.emMap.data() == 1))
    expected_at_2 = np.sum((map_levels.emMap.data() == 3)| (map_levels.emMap.data() == 2))
    expected_at_3 = np.sum(map_levels.emMap.data() == 3)
    assert vol_levels[1] == expected_at_3
    assert vol_levels[0.5] == expected_at_2
    assert vol_levels[0.3] == expected_at_1
    print('test_volume_isosurfaces OK')

# Test segment volumes at same isosurface
def test_volume_segments():
    map_segments = molecule.Molecule("maps/segments.map", 1)
    map_segments.generateSegments(3,1)
    segments_masks = map_segments.getSegmentsMasks()[1]
    vol_segments = map_segments.getSegmentsVolume()[1]
    for k in vol_segments.keys():
        assert vol_segments[k] == np.sum(segments_masks[k])
    print('test_volume_segments OK')

# Test correlation
def test_correlation_itself():
    map50 = molecule.Molecule('maps/cube50x50x50.map', 1)
    corr = map50.getCorrelation(map50)[1]
    assert corr == 1.0
    print('test_correlation_itself OK')

def test_correlation_bigger():
    map50 = molecule.Molecule('maps/cube50x50x50.map', 1)
    map100 = molecule.Molecule('maps/cube100x100x100.map', 1)
    corr = map50.getCorrelation(map100)[1]
    assert corr == 1.0
    print('test_correlation_bigger OK')

