"""
Unit tests for buildH_calcOP.

Test functions from module cli
"""

import pathlib
import pytest

import MDAnalysis as mda

from buildh import cli

DIR_DATA = "test_data"
PATH_ROOT_DATA = pathlib.Path(__file__).parent / DIR_DATA

# Ignore some MDAnalysis warnings for this test file
pytestmark = pytest.mark.filterwarnings('ignore::UserWarning')


class TestSlice():
    """
    Class for testing the slicing options from check_slice_options().

    This test trajectory contains:
       - 10 frames
       - from 0 to 10000 ps
    """

    # path for the Berger POPC files
    PATH_DATA = PATH_ROOT_DATA / "Berger_POPC"

    def setup_class(self):
        """Initialize attributes."""
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
        """Test with correct values.

        Parameters
        ----------
        begin : int or None
            the first frame to read (in ps). Can be 'None'.
        end : [type]
            the last frame to read (in ps). Can be 'None'.
        result : bool
            expected result
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
        Test if the correct exception is raised with wrong values.

        Parameters
        ----------
        begin : int or None
            the first frame to read (in ps).
        end : [type]
            the last frame to read (in ps).
        """
        with pytest.raises(IndexError) as err:
            cli.check_slice_options(self.universe, first_frame=begin, last_frame=end)
        assert "Incorrect slice options" in str(err.value)


class TestChecksAtoms():
    """Class for testing the 2 functions `check_def_file()` and `check_atom()`."""

    # path for the Berger POPC files
    PATH_DATA = PATH_ROOT_DATA / "Berger_POPC"

    def setup_class(self):
        """Initialize attributes."""
        pdb = self.PATH_DATA / "2POPC.pdb"

        self.universe = mda.Universe(str(pdb))


    @pytest.mark.parametrize('resname, name, result', [
        ("POPC","C1", True),
        ("POPC", "C48", True),
        ("POP", "C20", False),
        ("POPC", "C310", False)
    ])
    def test_check_atom(self, resname, name, result):
        """Test check_atom().

        Parameters
        ----------
        resname : str
            residue name of the atom
        name : str
            atom name to be tested.
        result : bool
            expected result
        """
        assert cli.check_atom(self.universe, resname, name) == result


    @pytest.mark.parametrize('resname, names, result, out', [
        ("POPC",["C1", "C38", "C5"], True, ""),
        ("POPC",["C1", "C38", "C310"], False, "Atom C310 of residue POPC"),
        ("POPC",["C211", "C38", "C5"], False, "Atom C211 of residue POPC"),
    ])
    def test_def_file(self, capfd, resname, names, result, out):
        """Test check_def_file().

        Parameters
        ----------
        capfd: pytest fixture
            Capture stdout/stderr output.
        resname : str
            residue name of the atoms
        names : list of str
            list of atom name to be tested.
        result : bool
            expected result
        out: str
            expected print output
        """
        res = cli.check_def_file(self.universe, resname, names)
        captured = capfd.readouterr()
        assert res == result
        assert out in captured.out
