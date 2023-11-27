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
        lower_bound = lower_bound[self.param_mask]
        upper_bound = upper_bound[self.param_mask]

        self.bsp = BSPtree(depth, lower_bound, upper_bound, idx)

    def idx(self, point):
        assert np.all(point[self.frozen_mask] == self.frozen_params)
        return self.bsp.at(point[self.param_mask])

    def at(self, point):
        assert np.all(point[self.frozen_mask] == self.frozen_params)
        return self.data[self.bsp.at(point[self.param_mask])]

    def sort(self, points):
        out = [[] for i in range(self.size)]
        for point in points:
            idx = self.idx(point)
            out[idx].append(point)

        return out
