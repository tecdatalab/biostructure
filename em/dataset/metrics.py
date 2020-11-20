import numpy as np


def intersection_over_union(segmented_map, gt_map):
    s_array = segmented_map.getEmMap().data()
    gt_array = gt_map.getEmMap().data()
    labels = np.unique(gt_array)
    if s_array.shape != gt_array.shape:
        return ValueError("Arrays must have same shape")
    else:
        iou_list = []
        for label in labels:
            if label == 0:
                continue
            s_mask = s_array == label
            gt_mask = gt_array == label
            overlap = np.sum(np.logical_and(s_mask,gt_mask))
            union = np.sum(np.logical_or(s_mask,gt_mask))
            iou_list.append(overlap/union)
        iou = np.mean(iou_list)
        return iou
            
def homogenity(segmented_map, gt_map):
    s_array = segmented_map.getEmMap().data()
    gt_array = gt_map.getEmMap().data()
    labels = np.unique(gt_array)
    if s_array.shape != gt_array.shape:
        return ValueError("Arrays must have same shape")
    else:
        h_list = []
        for label in labels:
            if label == 0:
                continue
            s_mask = s_array == label
            gt_mask = gt_array == label  
            overlap = np.sum(np.logical_and(s_mask,gt_mask))
            falses = np.sum(np.logical_xor(s_mask,gt_mask))
            h_list.append(falses/overlap)
        h = np.mean(h_list)
        return h

def proportion(segmented_map, gt_map):
    s_array = segmented_map.getEmMap().data()
    gt_array = gt_map.getEmMap().data()
    labels = np.unique(gt_array)
    if s_array.shape != gt_array.shape:
        return ValueError("Arrays must have same shape")
    else:
        p_list = []
        count_labels = 0
        for label in labels:
            if label == 0:
                continue
            s_mask = s_array == label
            gt_mask = gt_array == label
            # need to check if should remove 0s
            num_labels = len(np.unique(s_array[gt_mask]))
            p_list.append(num_labels)
            count_labels +=1
        p = np.sum(p_list)/count_labels
        return p
            
def consistency(segmented_map, gt_map):
    s_array = segmented_map.getEmMap().data()
    gt_array = gt_map.getEmMap().data()
    labels = np.unique(gt_array)
    if s_array.shape != gt_array.shape:
        return ValueError("Arrays must have same shape")
    else:
        c_list = []
        for label in labels:
            if label == 0:
                continue
            s_mask = s_array == label
            gt_mask = gt_array == label
            truth = np.sum(gt_mask)
            overlap = np.sum(np.logical_and(s_mask,gt_mask))
            c_list.append(overlap/truth)
        c = np.mean(c_list)
        return c 

