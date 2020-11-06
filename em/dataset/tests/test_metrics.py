from em.dataset.metrics import *
import em.molecule as molecule

def test_iou_same():
    segments_gt = molecule.Molecule('maps/segments_gt.map',1)
    iou = intersection_over_union(segments_gt, segments_gt)
    assert iou == 1.0

def test_iou_different():
    segments_gt = molecule.Molecule('maps/cubes_gt.map',1) 
    segments_bad_gt = molecule.Molecule('maps/cubes_gt_inv.map',1) 
    iou = intersection_over_union(segments_bad_gt, segments_gt)
    d1 = segments_gt.getEmMap().data()
    d2 = segments_bad_gt.getEmMap().data()
    labels1 = np.unique(d1)
    labels2 = np.unique(d2)
    labels_in_common = np.intersect1d(labels1, labels2)
    print(labels_in_common)
    # Remove label 0
    labels_in_common = np.delete(labels_in_common, np.argwhere(labels_in_common == 0))
    expected_iou = 0
    for label in labels_in_common:
        intersec = np.sum(np.logical_and(d1==label,d2==label))
        union = np.sum(np.logical_or(d1==label,d2==label))
        expected_iou +=intersec/union
    expected_iou /= len(labels_in_common) 
    assert iou == expected_iou
