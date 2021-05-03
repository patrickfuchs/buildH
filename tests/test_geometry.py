"""Unit tests for buildH_calcOP.

Test functions from module geometry.

Notes
-----
All these values were extracted on a simple example test_data/10POPC.pdb.
Beware, we realized that sometimes print() displays different values than the
real ones, which could make the tests fail when using assert_almost_equal()
(probably a rounding problem). To get the right value, use instead a debugger
(e.g. in VScode) or test the function in the Python interpreter.
"""

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

from buildh import geometry as geom

@pytest.mark.parametrize('vec, result', [
    (np.array([-1.165698, 1.3688029, -0.6189914]), 1.9014794),
])
def test_norm(vec, result):
    """Test norm().

    Parameters
    ----------
    vec : numpy 1D-array
        input data
    result : float
        reference result
    """
    assert_almost_equal(geom.norm(vec), result)


@pytest.mark.parametrize('vec, result', [
    (np.array([-1.165698, 1.3688029, -0.6189914]),
     np.array([-0.61304796, 0.7198621, -0.32553148])
     ),
])
def test_normalize(vec, result):
    """Test normalize().

    Parameters
    ----------
    vec : numpy 1D-array
        input data
    result : numpy 1D-array
        reference result
    """
    assert_almost_equal(geom.normalize(vec), result)


@pytest.mark.parametrize('atom1, atom2, atom3, result', [
    (np.array([21.13, 41.14, 31.36]), np.array([20.64, 39.8, 31.91]),
     np.array([20.4, 39.52, 33.25]), 2.18789998),
])
def test_calc_angle(atom1, atom2, atom3, result):
    """Test calc_angle().

    Parameters
    ----------
    atom1 : numpy 1D-array
        coordinates
    atom2 : numpy 1D-array
        coordinates
    atom3 : numpy 1D-array
        coordinates
    result : float
        reference angle
    """
    assert_almost_equal(geom.calc_angle(atom1, atom2, atom3), result)


@pytest.mark.parametrize('vec, theta, result', [
    (np.array([-0.61304796, 0.7198621, -0.32553148]), 1.9106332362490186,
     np.array([0.57735027, -0.50055158, 0.58776498, -0.26579535]),
     ),
])
def test_vec2quaternion(vec, theta, result):
    """Test vec2quaternion().

    Parameters
    ----------
    vec : numpy 1D-array
        input data
    theta: float
        input angle
    result : numpy 1D-array
        reference quaternion
    """
    assert_almost_equal(geom.vec2quaternion(vec, theta), result)


@pytest.mark.parametrize('quaternion, result', [
    (np.array([0.57735027, -0.50055158, 0.58776498, -0.26579535]),
     np.array([[ 0.16777038, -0.28149935,  0.9447811],
               [-0.89532741,  0.35760195,  0.26553678],
               [-0.41260397, -0.89043758, -0.19203905]]),
     ),
])
def test_calc_rotation_matrix(quaternion, result):
    """Test calc_rotation_matrix().

    Parameters
    ----------
    quaternion : numpy 1D-array
        input data
    result : numpy 2D-array
        reference matrix
    """
    assert_almost_equal(geom.calc_rotation_matrix(quaternion), result)


@pytest.mark.parametrize('vec_to_rotate, rotation_axis, rad_angle, result', [
    (np.array([-1.0199966, -0.4300003,  0.9700012], dtype=np.float32),
     np.array([-0.61304796,  0.7198621,  -0.32553148]),
     1.9106332362490186,
     np.array([0.58863854, 0.69101293, 0.41953045])
     ),
])
def test_apply_rotation(vec_to_rotate, rotation_axis, rad_angle, result):
    """Test apply_rotation().

    Parameters
    ----------
    vec_to_rotate : numpy 1D-array
        input data
    rotation_axis : numpy 1D-array
        input data
    rad_angle: float
        input angle
    result : numpy 1D-array
        reference matrix
    """
    assert_almost_equal(geom.apply_rotation(vec_to_rotate, rotation_axis, rad_angle), result)


@pytest.mark.parametrize('A, B, result', [
    (np.array([-2.0, -1.449997, 0.5600014]),
     np.array([-1.0199966, -0.4300003, 0.9700012]),
     np.array([-1.165698, 1.3688029, -0.6189914])
     ),
])
def test_cross_product(A, B, result):
    """Test cross_product().

    Parameters
    ----------
    A : numpy 1D-array
        input data
    B : numpy 1D-array
        input data
    result : numpy 1D-array
        reference result
    """
    assert_almost_equal(geom.cross_product(A, B), result)


@pytest.mark.parametrize('C, H, result', [
    (np.array([34.42, 46.94, 26.31]),
     np.array([35.06161421, 47.69320272, 26.76728762]),
    -0.23599087203193325
     ),
])
def test_calc_OP(C, H, result):
    """Test calc_OP().

    Parameters
    ----------
    C : numpy 1D-array
        input data
    H : numpy 1D-array
        input data
    result : float
        reference angle
    """
    assert_almost_equal(geom.calc_OP(C, H), result)
