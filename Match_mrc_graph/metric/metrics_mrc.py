import numpy as np


def get_geometric_overlap_p(mrc_data1, mrc_data2):
    mrc1_points_tuple = np.where(mrc_data1 > 0)
    mrc2_points_tuple = np.where(mrc_data2 > 0)
    mass1 = mrc1_points_tuple[0].shape[0]
    mass2 = mrc2_points_tuple[0].shape[0]
    tuples_mrc1 = np.column_stack(mrc1_points_tuple)
    tuples_mrc2 = np.column_stack(mrc2_points_tuple)

    all_points = np.concatenate((tuples_mrc1, tuples_mrc2), axis=0)
    unique_points = np.unique(all_points, axis=0)  # This is the slow part
    mass_overlap = all_points.shape[0] - unique_points.shape[0]
    result = mass_overlap / min(mass1, mass2)
    return mass1, mass2, mass_overlap, result


def get_cross_correlation(mrc_data1, mrc_data2):
    CC = np.corrcoef(mrc_data1.ravel(), mrc_data2.ravel())
    result = CC[0, 1]
    return result
