import pdb

import scipy as sp
import warnings

def filterBid(allids, sbids):
    '''
    gets two id list
    returns matching index
    '''

    if sp.unique(sbids).shape[0]  != sbids.shape[0]:
        warnings.warn("superset ids are not unique: Making it unique")
        sbids = sp.unique(sbids)
    if sp.unique(allids).shape[0] != allids.shape[0]:
        warnings.warn("Subset ids are not unique: Making it unique")
        allids = sp.unique(allids)
    if sp.sum(sp.sort(allids) == allids) != allids.shape[0]: 
        warnings.warn("Superset ids are not sorted: Sorting it")
        allids = sp.sort(allids)
    if sp.sum(sp.sort(sbids) == sbids) != sbids.shape[0]:
        warnings.warn('subset ids are not sorted: Sorting it')
        sbids = sp.sort(sbids)
    return sp.where(sp.in1d(allids, sbids))[0]
  


def unique_rows_idx(a, return_counts = False):

    unique_a = sp.unique(a.view([('', a.dtype)]*a.shape[1]),return_index = True, return_counts = return_counts)
    if return_counts:
        return unique_a[1], unique_a[2]
    return unique_a[1]



def unique_rows(array, index = None, counts = False):
    """Make array unique by rows"""

    if array.shape[0] == 0:
        if index == True:
            return (array, [])
        else:
            return (array)

    if len(array.shape) == 1:
        if index == True:
            return (array, [0])
        else:
            return array

    (array_s, s_idx) = sort_rows(array, True)
    tmp = [False]
    tmp.extend([sp.all(array_s[i-1, :] == array_s[i, :]) for i in range(1, array.shape[0])])
    k_idx = sp.where(~sp.array(tmp, dtype='bool'))[0]

    if index == True:
        if counts:
            uidx = s_idx[k_idx]
            dist = sp.append((uidx[1:] - uidx[:-1]), array.shape[0] - uidx[-1])
            return (array[s_idx[k_idx], :], s_idx[k_idx], dist)
        else:
            return (array[s_idx[k_idx], :], s_idx[k_idx])
    else:
        return array[s_idx[k_idx], :]

def sort_rows(array, index = None):
    """Sort array by rows"""

    ### empty array
    if array.shape[0] == 0:
        if index == True:
            return (array, [])
        else:
            return (array)

    ### only one row
    if len(array.shape) == 1:
        if index == True:
            return (array, [0])
        else:
            return (array)

    ### more than one row
    s_idx = sp.lexsort([array[:, -i] for i in range(1, array.shape[1] + 1)])

    if index == True:
        return (array[s_idx, :], s_idx)
    else:
        return array[s_idx, :]

