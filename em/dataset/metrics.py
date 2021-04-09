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
            #if label == 0:
            #    continue
            s_mask = s_array == label
            gt_mask = gt_array == label
            overlap = np.sum(np.logical_and(s_mask,gt_mask))
            union = np.sum(np.logical_or(s_mask,gt_mask))
            iou_list.append(overlap/union)
        iou = np.mean(iou_list)
        return iou


def average_precision(segmented_map, gt_map, thresholds=np.arange(0.5,1,0.05)):
    segmented_array = segmented_map.getEmMap().data()
    gt_array = gt_map.getEmMap().data()
    
    segmented_labels = np.unique(segmented_array)
    segmented_labels = segmented_labels[segmented_labels!=0]
    gt_labels = np.unique(gt_array)
    gt_labels = gt_labels[gt_labels!=0]
    segmented_masks = [ segmented_array==l for l in segmented_labels ]
    gt_masks = [ gt_array==l for l in gt_labels ]
   
    iou_tensor = np.zeros([len(thresholds), len(segmented_masks), len(gt_masks)])
    for i,seg_mask in enumerate(segmented_masks):
        for j,gt_mask in enumerate(gt_masks):
            iou_tensor[:, i, j] = iou_at_thresholds(gt_mask, seg_mask, thresholds)
    TP = np.sum((np.sum(iou_tensor, axis=2) == 1), axis=1)
    FP = np.sum((np.sum(iou_tensor, axis=1) == 0), axis=1)
    FN = np.sum((np.sum(iou_tensor, axis=2) == 0), axis=1)

    precision = TP / (TP + FP + FN)

    return np.mean(precision)

def iou_at_thresholds(gt_mask, seg_mask, thresholds=np.arange(0.5,1,0.05)):
    intersection = np.logical_and(gt_mask, seg_mask)
    union = np.logical_or(gt_mask, seg_mask)
    iou = np.sum(intersection > 0) / np.sum(union > 0)
    return iou > thresholds


def dice(segmented_map, gt_map):
    s_array = segmented_map.getEmMap().data()
    gt_array = gt_map.getEmMap().data()
    labels = np.unique(gt_array)
    if s_array.shape != gt_array.shape:
        return ValueError("Arrays must have same shape")
    else:
        dice_list = []
        for label in labels:
            #if label == 0:
            #    continue
            s_mask = s_array == label
            gt_mask = gt_array == label
            overlap = np.sum(np.logical_and(s_mask,gt_mask))
            added = np.sum(s_mask) + np.sum(gt_mask)
            dice_list.append(2*overlap/added)
        dice = np.mean(dice_list)
        return dice
            
def homogenity(segmented_map, gt_map):
    s_array = segmented_map.getEmMap().data()
    gt_array = gt_map.getEmMap().data()
    labels = np.unique(s_array)
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
            volume = np.sum(gt_mask)
            h_list.append(overlap/(volume+np.finfo(float).eps))
            print("label {} overlap {} falses {} result {}".format(label, overlap,volume, overlap/(volume+np.finfo(float).eps)))
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
            s_mask_non_background = s_array != 0
            proportion_mask = gt_mask * s_mask_non_background
            # need to check if should remove 0s
            num_labels = len(np.unique(s_array[proportion_mask]))
            p_list.append(num_labels)
            count_labels +=1
            print("label {} proportion {}".format(label, num_labels))
        p = np.sum(p_list)/count_labels
        return p
            
def consistency(segmented_map, gt_map):
    s_array = segmented_map.getEmMap().data()
    gt_array = gt_map.getEmMap().data()
    labels = np.unique(gt_array)
    if s_array.shape != gt_array.shape:
        return ValueError("Arrays must have same shape")
    else:
        volumes_dict = {}
        for label in labels:
            if label == 0:
                continue
            s_mask = s_array == label
            gt_mask = gt_array == label
            volumes_dict[label] = np.sum(s_mask)
        label = max(volumes_dict, key=volumes_dict.get)
        gt_mask = gt_array == label
        c = len(np.unique(s_array[gt_mask])) 
        return c 

