#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""
Unit tests for buildH_calcOP

Test functions from module hydrogens
"""

import collections
import numpy as np
from numpy.testing import (
    assert_almost_equal,
)
import pytest

import hydrogens


class TestComputeHydrogen:
    """
    Test class for the functions which reconstruct hydrogens from atom and helpers
    """

    # Data tuple for get_CH() test
    Data = collections.namedtuple('Data', ['atom', 'helper1', 'helper2', 'helper3', 'H_coord'])

    @pytest.mark.parametrize('data',[
        (Data(np.array([26.87, 45.09, 26.03], dtype=np.float32),
              np.array([27.85, 46.07, 25.38], dtype=np.float32),
              np.array([27.65, 44.79, 27.32], dtype=np.float32),
              np.array([25.53, 45.59, 26.01], dtype=np.float32),
              np.array([26.61868981, 44.14296091, 25.552445]))),
        (Data(np.array([09.37, 46.61, 21.07], dtype=np.float32),
              np.array([08.62, 47.02, 19.80], dtype=np.float32),
              np.array([09.71, 47.93, 21.77], dtype=np.float32),
              np.array([08.64, 45.68, 21.88], dtype=np.float32),
              np.array([10.27737988, 46.049353, 20.84542051]))),
        (Data(np.array([30.83, 37.09, 27.29], dtype=np.float32),
              np.array([30.72, 36.37, 25.95], dtype=np.float32),
              np.array([31.63, 38.33, 26.88], dtype=np.float32),
              np.array([29.69, 37.27, 28.14], dtype=np.float32),
              np.array([31.30033387, 36.45885357, 28.04401681]))),
        (Data(np.array([40.79, 32.48, 25.71], dtype=np.float32),
              np.array([41.58, 32.30, 24.41], dtype=np.float32),
              np.array([39.48, 31.71, 25.54], dtype=np.float32),
              np.array([41.50, 31.96, 26.85], dtype=np.float32),
              np.array([40.62567756, 33.54206805, 25.89195599]))),
    ])
    def test_get_CH(self, data):
        """
        Test for get_CH()
        3 helpers are needed and 1 hydrogen is returned
        """
        assert_almost_equal(hydrogens.get_CH(data.atom, data.helper1, data.helper2, data.helper3),
                            data.H_coord)


    # Data tuple for get_CH2() test
    Data = collections.namedtuple('Data', ['atom', 'helper1', 'helper2', 'H1_coord', 'H2_coord'])

    @pytest.mark.parametrize('data',[
        (Data(np.array([32.42, 45.49, 26.87], dtype=np.float32),
              np.array([33.40, 46.51, 27.28], dtype=np.float32),
              np.array([31.74, 45.57, 25.50], dtype=np.float32),
              np.array([31.62418725, 45.5047426 , 27.61469385]),
              np.array([32.93474471, 44.52992542, 26.90727792]))),
        (Data(np.array([28.81, 46.72, 30.22], dtype=np.float32),
              np.array([28.65, 45.64, 29.15], dtype=np.float32),
              np.array([28.1 , 46.33, 31.52], dtype=np.float32),
              np.array([28.38370523, 47.65215101, 29.84922984]),
              np.array([29.87120526, 46.8618091 , 30.42452995]))),
        (Data(np.array([78.72, 55.93, 30.95], dtype=np.float32),
              np.array([77.34, 55.74, 31.57], dtype=np.float32),
              np.array([78.85, 55.49, 29.49], dtype=np.float32),
              np.array([79.43524658, 55.35603092, 31.53913845]),
              np.array([78.97064378, 56.98933657, 31.0055434 ]))),
        (Data(np.array([15.51, 14.20, 40.00], dtype=np.float32),
              np.array([14.35, 14.44, 39.03], dtype=np.float32),
              np.array([15.67, 15.44, 40.89], dtype=np.float32),
              np.array([16.42873323, 14.03015964, 39.43858751]),
              np.array([15.29632013, 13.32883868, 40.61928919]))),
    ])
    def test_get_CH2(self, data):
        """
        Test for get_CH2()
        2 helpers are needed and 2 hydrogens are returned
        """
        assert_almost_equal(hydrogens.get_CH2(data.atom, data.helper1, data.helper2),
                            (data.H1_coord, data.H2_coord))


    # Data tuple for get_CH3() test
    Data = collections.namedtuple('Data', ['atom', 'helper1', 'helper2', 'H1_coord', 'H2_coord', 'H3_coord'])

    @pytest.mark.parametrize('data',[
        (Data(np.array([34.42, 46.94, 26.31], dtype=np.float32),
              np.array([33.40, 46.51, 27.28], dtype=np.float32),
              np.array([32.42, 45.49, 26.87], dtype=np.float32),
              np.array([35.06161421, 47.69320272, 26.76728762]),
              np.array([33.93128850, 47.36328732, 25.43245201]),
              np.array([35.02249090, 46.08195972, 26.01188584]))),
        (Data(np.array([31.05, 43.87, 46.61], dtype=np.float32),
              np.array([30.47, 43.67, 45.21], dtype=np.float32),
              np.array([30.13, 45.04, 44.62], dtype=np.float32),
              np.array([31.29707701, 42.90060506, 47.04281476]),
              np.array([30.31528839, 44.37225161, 47.23931947]),
              np.array([31.95123226, 44.47976139, 46.54621353]))),
        (Data(np.array([63.69, 53.48, 16.38], dtype=np.float32),
              np.array([64.97, 52.74, 16.35], dtype=np.float32),
              np.array([65.30, 52.00, 17.57], dtype=np.float32),
              np.array([63.54372877, 53.99164339, 15.42872333]),
              np.array([63.71263550, 54.21253590, 17.18683127]),
              np.array([62.87017354, 52.78125334, 16.54655418]))),
        (Data(np.array([05.02, 23.08, 20.84], dtype=np.float32),
              np.array([04.46, 24.37, 21.26], dtype=np.float32),
              np.array([03.10, 24.66, 20.78], dtype=np.float32),
              np.array([06.02425455, 22.96939049, 21.2490702 ]),
              np.array([05.06460696, 23.04033616, 19.75163578]),
              np.array([04.38703115, 22.27223366, 21.20737478]))),
    ])
    def test_get_CH3(self, data):
        """
        Test for get_CH3()
        2 helpers are needed and 3 hydrogens are returned
        """
        assert_almost_equal(hydrogens.get_CH3(data.atom, data.helper1, data.helper2),
                            (data.H1_coord, data.H2_coord, data.H3_coord))


    # Data tuple for get_CH_double_bond() test
    Data = collections.namedtuple('Data', ['atom', 'helper1', 'helper2', 'H1_coord'])

    @pytest.mark.parametrize('data',[
        (Data(np.array([20.64, 39.80, 31.91], dtype=np.float32),
              np.array([21.13, 41.14, 31.36], dtype=np.float32),
              np.array([20.40, 39.52, 33.25], dtype=np.float32),
              np.array([20.46454439, 38.99866187, 31.19224364]))),
        (Data(np.array([15.41, 13.27, 34.47], dtype=np.float32),
              np.array([15.05, 14.68, 34.00], dtype=np.float32),
              np.array([15.71, 12.86, 35.77], dtype=np.float32),
              np.array([15.43518964, 12.49687463, 33.70205465]))),
        (Data(np.array([07.29, 24.63, 35.47], dtype=np.float32),
              np.array([07.65, 24.12, 34.08], dtype=np.float32),
              np.array([07.81, 25.78, 36.05], dtype=np.float32),
              np.array([06.57262835, 24.05054190, 36.05112706]))),

    ])
    def test_get_CH_double_bond(self, data):
        """
        Test for get_CH_double_bond()
        2 helpers are needed and 1 hydrogen is returned
        """
        assert_almost_equal(hydrogens.get_CH_double_bond(data.atom, data.helper1, data.helper2),
                            data.H1_coord)
