"""
Unit tests for buildH_calcOP

Test functions from module cli
"""

import pathlib
import pytest

import MDAnalysis as mda

from buildh import cli

dir_data = "test_data"
path_root_data = pathlib.Path(__file__).parent / dir_data

# Ignore some MDAnalysis warnings for this test file
pytestmark = pytest.mark.filterwarnings('ignore::UserWarning')


class TestSlice():
    """
    Class for testing the slicing options from check_slice_options()
    This test trajectory contains:
       - 10 frames
       - from 0 to 10000 ps
    """
    # path for the Berger POPC files
    PATH_DATA = path_root_data / "Berger_POPC"

    def setup_class(self):
        pdb = self.PATH_DATA / "2POPC.pdb"
        xtc = self.PATH_DATA / "2POPC.xtc"

        self.universe = mda.Universe(str(pdb), str(xtc))

    # 'begin' and 'end' are in ps.
    # result_ref is a tuple containing the corresponding frame number for the slice option.
    @pytest.mark.parametrize('begin, end, result', [
        (None, None,  (0,11)),
        (None, 3000,  (0, 4)),
        (4000, None,  (4, 11)),
        (6000, 10000, (6,11)),
        (0,    9501,  (0, 11)),
        (0,    9499,  (0, 10)),
        (2501, 3000,  (3,4)),
        (2499, 3000,  (2,4)),
        (7501, 7550,  (8,9)),
    ])
    def test_correct_slice(self, begin, end, result):
        """
        Test with correct values
        """
        assert cli.check_slice_options(self.universe, first_frame=begin, last_frame=end) == result

    @pytest.mark.parametrize('begin, end', [
        (0,     -1000),
        (-1000, 2000),
        (400,   20000),
        (2000,  1000),
        (11000, 12000),
    ])
    def test_wrong_slice(self, begin, end):
        """
        Test if the correct exception is raised
        """

        with pytest.raises(IndexError) as err:
             cli.check_slice_options(self.universe, first_frame=begin, last_frame=end)
        assert "Incorrect slice options" in str(err.value)

