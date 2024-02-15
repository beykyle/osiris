import osiris
from unittest import TestCase
import numpy as np


class BSPTest(TestCase):
    def test_bsp_int(self):
        bsp = osiris.BinarySPTree(
            2,
            np.array([-1, -1, -1]),
            np.array([1, 1, 1]),
            np.array([0, 1, 2, 3]),
        )

        np.testing.assert_equal(bsp.at(np.array([-1.0 / 2, -1.0 / 2, -1.0 / 2])), 0)
        np.testing.assert_equal(bsp.at(np.array([-1.0 / 2, -1.0 / 2, 1.0 / 2])), 0)

        np.testing.assert_equal(bsp.at(np.array([-1.0 / 2, 1.0 / 2, -1.0 / 2])), 1)
        np.testing.assert_equal(bsp.at(np.array([-1.0 / 2, 1.0 / 2, 1.0 / 2])), 1)

        np.testing.assert_equal(bsp.at(np.array([1.0 / 2, -1.0 / 2, -1.0 / 2])), 2)
        np.testing.assert_equal(bsp.at(np.array([1.0 / 2, -1.0 / 2, 1.0 / 2])), 2)

        np.testing.assert_equal(bsp.at(np.array([1.0 / 2, 1.0 / 2, -1.0 / 2])), 3)
        np.testing.assert_equal(bsp.at(np.array([1.0 / 2, 1.0 / 2, 1.0 / 2])), 3)

    def test_bsp_real(self):
        bsp = osiris.BinarySPTree(
            2,
            np.array([-1, -1, -1]),
            np.array([1, 1, 1]),
            np.array([0.0, 0.1, 0.2, 0.3], dtype=object),
        )

        np.testing.assert_equal(bsp.at(np.array([-1.0 / 2, -1.0 / 2, -1.0 / 2])), 0)
        np.testing.assert_equal(bsp.at(np.array([-1.0 / 2, -1.0 / 2, 1.0 / 2])), 0)

        np.testing.assert_equal(bsp.at(np.array([-1.0 / 2, 1.0 / 2, -1.0 / 2])), 0.1)
        np.testing.assert_equal(bsp.at(np.array([-1.0 / 2, 1.0 / 2, 1.0 / 2])), 0.1)

        np.testing.assert_equal(bsp.at(np.array([1.0 / 2, -1.0 / 2, -1.0 / 2])), 0.2)
        np.testing.assert_equal(bsp.at(np.array([1.0 / 2, -1.0 / 2, 1.0 / 2])), 0.2)

        np.testing.assert_equal(bsp.at(np.array([1.0 / 2, 1.0 / 2, -1.0 / 2])), 0.3)
        np.testing.assert_equal(bsp.at(np.array([1.0 / 2, 1.0 / 2, 1.0 / 2])), 0.3)

    def test_bsp_string(self):
        bsp = osiris.BinarySPTree(
            2,
            np.array([-1, -1, -1]),
            np.array([1, 1, 1]),
            np.array(["a", "b", "c", "d"], dtype=object),
        )

        np.testing.assert_string_equal(
            bsp.at(np.array([-1.0 / 2, -1.0 / 2, -1.0 / 2])), "a"
        )
        np.testing.assert_string_equal(
            bsp.at(np.array([-1.0 / 2, -1.0 / 2, 1.0 / 2])), "a"
        )

        np.testing.assert_string_equal(
            bsp.at(np.array([-1.0 / 2, 1.0 / 2, -1.0 / 2])), "b"
        )
        np.testing.assert_string_equal(
            bsp.at(np.array([-1.0 / 2, 1.0 / 2, 1.0 / 2])), "b"
        )

        np.testing.assert_string_equal(
            bsp.at(np.array([1.0 / 2, -1.0 / 2, -1.0 / 2])), "c"
        )
        np.testing.assert_string_equal(
            bsp.at(np.array([1.0 / 2, -1.0 / 2, 1.0 / 2])), "c"
        )

        np.testing.assert_string_equal(
            bsp.at(np.array([1.0 / 2, 1.0 / 2, -1.0 / 2])), "d"
        )
        np.testing.assert_string_equal(
            bsp.at(np.array([1.0 / 2, 1.0 / 2, 1.0 / 2])), "d"
        )

    def test_bsp_equal_bounds(self):
        # This bsp is only 2D; the last dimension will be ignored bc the upper and lower
        # bounds are equal.
        # This means we split x at 0, then y at 0, then x at +/- 1/2
        # divides 2D box into 8 sub-boxes using planes y=0 and x=0, x=-1/2, x=+1/2
        #
        #     y
        #     |
        #   1 -------------
        #     |2 |3 |6 |7 |
        #   0 -------------
        #     |0 |1 |4 |5 |
        #  -1 ---------------x
        #    -1 -1/2  1/2  1
        #
        bsp = osiris.BinarySPTree(
            3,
            np.array([-1, -1, 1]),
            np.array([1, 1, 1]),
            np.arange(0, 8, 1),
        )

        np.testing.assert_equal(bsp.at(np.array([-0.75, -0.75, 1])), 0)
        np.testing.assert_equal(bsp.at(np.array([-0.75, -0.75, 1])), 0)

        np.testing.assert_equal(bsp.at(np.array([0.75, 0.75, 1])), 7)
        np.testing.assert_equal(bsp.at(np.array([0.75, 0.75, 1])), 7)

    def test_sort(self):
        bsp = osiris.BinarySPTree(
            2,
            np.array([-1, -1, -1]),
            np.array([1, 1, 1]),
            np.array([0, 1, 2, 3]),
        )

        params = np.array([[-0.5, -0.5, -0.5], [0.9, 0.9, 0.9], [0.8, 0.8, 0.8]])
        sorted_params = bsp.sort(params)

        np.testing.assert_equal(sorted_params[0], [0])
        np.testing.assert_equal(sorted_params[1], [])
        np.testing.assert_equal(sorted_params[2], [])
        np.testing.assert_equal(sorted_params[3], [1,2])

    def test_sort_mask(self):
        bsp = osiris.BinarySPTree(
            2,
            np.array([-1, -1, -1, 8]),
            np.array([1, 1, 1, 8]),
            np.array([0, 1, 2, 3]),
        )

        params = np.array(
            [[1, -0.5, -0.5, -0.5, 8], [8, 0.9, 0.9, 0.9, 8], [9, 0.8, 0.8, 0.8, 8]]
        )
        sorted_params = bsp.sort(params, mask=np.array([0, 1, 1, 1, 1], dtype=bool))

        np.testing.assert_equal(sorted_params[0], [0])
        np.testing.assert_equal(sorted_params[1], [])
        np.testing.assert_equal(sorted_params[2], [])
        np.testing.assert_equal(sorted_params[3], [1,2])

    def test_sp_bounds(self):
        #     y
        #     |
        #   1 -------
        #     |1 |3 |
        #   0 -------
        #     |0 |2 |
        #  -1 --------x
        #    -1  0   1
        #
        bsp = osiris.BinarySPTree(
            2,
            np.array([-1, -1, -1]),
            np.array([1, 1, 1]),
            np.array([0, 1, 2, 3]),
        )

        expected_lbounds = np.array(
            [[-1, -1, -1], [-1, 0, -1], [0, -1, -1], [0, 0, -1]]
        )
        expected_ubounds = np.array([[0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1]])

        lbounds, ubounds = bsp.get_sub_partition_bounds()

        np.testing.assert_equal(lbounds, expected_lbounds)
        np.testing.assert_equal(ubounds, expected_ubounds)
