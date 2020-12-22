"""
Unit tests for buildH_calcOP

Test functions from module geomtry
"""

import numpy as np
from numpy.testing import (
    assert_almost_equal,
)
import pytest

from buildh import geometry


@pytest.mark.parametrize('atom1, atom2, atom3, result', [
    (np.array([21.13, 41.14, 31.36]), np.array([20.64, 39.8, 31.91]),
     np.array([20.4, 39.52, 33.25]), 2.18789998),
])
def test_angle(atom1, atom2, atom3, result):
    assert_almost_equal(geometry.calc_angle(atom1, atom2, atom3), result)
