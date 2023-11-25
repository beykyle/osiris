from .osiris_core import BinarySPTree_int as BSPtree
import numpy as np


class BinarySPTree:
    def __init__(self, depth, lower_bound, upper_bound, data):
        self.data = data
        idx = np.arange(0, len(data))
        self.bsp = BSPtree(depth, lower_bound, upper_bound, idx)

    def idx(self, point):
        return self.bsp.at(point)

    def at(self, point):
        return self.data[self.bsp.at(point)]
