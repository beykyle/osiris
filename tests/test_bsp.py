import osiris
from unittest import TestCase
import numpy as np


class BSPTest(TestCase):
    def test_bsp_simple(self):
        bsp = osiris.BinarySPTree(
            2, np.array([-1, -1, -1]), np.array([1, 1, 1]), np.array([0, 1, 2, 3])
        )

        np.testing.assert_equal(bsp(np.array([-1.0 / 2, -1.0 / 2, -1.0 / 2])), 0)
        np.testing.assert_equal(bsp(np.array([-1.0 / 2, -1.0 / 2, 1.0 / 2])), 0)

        np.testing.assert_equal(bsp(np.array([-1.0 / 2, 1.0 / 2, -1.0 / 2])), 1)
        np.testing.assert_equal(bsp(np.array([-1.0 / 2, 1.0 / 2, 1.0 / 2])), 1)

        np.testing.assert_equal(bsp(np.array([1.0 / 2, -1.0 / 2, -1.0 / 2])), 2)
        np.testing.assert_equal(bsp(np.array([1.0 / 2, -1.0 / 2, 1.0 / 2])), 2)

        np.testing.assert_equal(bsp(np.array([1.0 / 2, 1.0 / 2, -1.0 / 2])), 3)
        np.testing.assert_equal(bsp(np.array([1.0 / 2, 1.0 / 2, 1.0 / 2])), 3)
