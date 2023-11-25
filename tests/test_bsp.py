import osiris
from unittest import TestCase
import numpy as np


class BSPTest(TestCase):
    def test_bsp_int(self):
        bsp = osiris.BinarySPTree_int(
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

    def test_bsp_real(self):
        bsp = osiris.BinarySPTree_double(
            2,
            np.array([-1, -1, -1]),
            np.array([1, 1, 1]),
            np.array([0.0, 0.1, 0.2, 0.3], dtype=object),
        )

        np.testing.assert_equal(bsp(np.array([-1.0 / 2, -1.0 / 2, -1.0 / 2])), 0)
        np.testing.assert_equal(bsp(np.array([-1.0 / 2, -1.0 / 2, 1.0 / 2])), 0)

        np.testing.assert_equal(bsp(np.array([-1.0 / 2, 1.0 / 2, -1.0 / 2])), 0.1)
        np.testing.assert_equal(bsp(np.array([-1.0 / 2, 1.0 / 2, 1.0 / 2])), 0.1)

        np.testing.assert_equal(bsp(np.array([1.0 / 2, -1.0 / 2, -1.0 / 2])), 0.2)
        np.testing.assert_equal(bsp(np.array([1.0 / 2, -1.0 / 2, 1.0 / 2])), 0.2)

        np.testing.assert_equal(bsp(np.array([1.0 / 2, 1.0 / 2, -1.0 / 2])), 0.3)
        np.testing.assert_equal(bsp(np.array([1.0 / 2, 1.0 / 2, 1.0 / 2])), 0.3)
