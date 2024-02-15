from .osiris_core import BinarySPTree_int as BSPtree
import numpy as np


class BinarySPTree:
    def __init__(self, depth, lower_bound, upper_bound, data):
        self.depth = depth
        self.size = 2**depth
        self.data = data
        self.dimensions = lower_bound.size
        idx = np.arange(0, len(data))

        # first sanitizing bounds to freeze parameters that are equal
        mask = lower_bound == upper_bound
        self.frozen_param_idxs = np.nonzero(mask)
        self.frozen_params = lower_bound[mask]
        n_unfrozen = lower_bound[np.logical_not(mask)].size
        self.active_dimensions = n_unfrozen

        # bounds for unfrozen params only
        self.param_mask = np.logical_not(mask)
        self.frozen_mask = mask
        self.lower_bound = lower_bound[self.param_mask]
        self.upper_bound = upper_bound[self.param_mask]

        self.bsp = BSPtree(depth, self.lower_bound, self.upper_bound, idx)

    def idx(self, point):
        assert np.all(point[self.frozen_mask] == self.frozen_params)
        return self.bsp.at(point[self.param_mask])

    def at(self, point):
        assert np.all(point[self.frozen_mask] == self.frozen_params)
        return self.data[self.bsp.at(point[self.param_mask])]

    def sort(self, points, mask=None):
        """
        Params:
            points (ndarray): 2d array (mxn), m being number of points, n being the dimension
            mask (ndarray): logical mask of size n (the dimension of the points), enabling sorted a
                set of points that are partitioned alng only a subset of their dimensions
        Returns:
            out (list): list with an element for each partition containing the indices
                (along axis 0 of points) corresponding to the points w/in the partition
        """
        if mask is None:
            mask = np.ones(self.dimensions, dtype=np.int32)
        assert np.sum(mask) == self.dimensions
        out = [[] for i in range(self.size)]
        for i in range(points.shape[0]):
            p = points[i, mask]
            idx = self.idx(p)
            out[idx].append(i)
        out = [np.array(sv) for sv in out]
        return out

    def get_sub_partition_bounds(self):
        r"""
            Returns :
        2 arrays indxed by (number of sub partitions, number of dimensions),
        the first representing the lower bound of each partition, and the
        second the upper
        """
        lower_bounds = np.empty((self.size, self.dimensions))
        upper_bounds = np.empty((self.size, self.dimensions))
        for i in range(self.size):
            l, u = self.bsp.get_bounds(i)
            lower_bounds[i, self.param_mask] = l
            upper_bounds[i, self.param_mask] = u
            lower_bounds[i, self.frozen_mask] = self.frozen_params
            upper_bounds[i, self.frozen_mask] = self.frozen_params

        return lower_bounds, upper_bounds
