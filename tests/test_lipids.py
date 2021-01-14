"""
Unit tests for buildH_calcOP.

Test functions from module lipids
"""

import pathlib
import pytest

import MDAnalysis as mda

from buildh import lipids

dir_data = "test_data"
path_data = pathlib.Path(__file__).parent / dir_data

# Ignore some MDAnalysis warnings for this test file
pytestmark = pytest.mark.filterwarnings('ignore::UserWarning')


class TestLipids():
    """Test class for functions in the module lipids.py."""

    def setup_class(self):
        """Initialize attributes."""
        filenames = ["Berger_POPC.json", "CHARMM_POPC.json"]
        self.path_files = [path_data / f for f in filenames]

        pdb = path_data / "2POPC.pdb"
        self.universe = mda.Universe(str(pdb))


    def test_check_read_lipids_topH_success(self):
        """Test check_read_lipids_topH() with correct files."""
        lipids_tops = lipids.read_lipids_topH(self.path_files)

        # key : key of the outer dict
        # values : number of atoms in each dict
        reference = {'Berger_POPC': 41, 'Berger_PLA': 41, 'Berger_POP':41, 'CHARMM_POPC':41}

        assert len(lipids_tops) == 4
        assert lipids_tops.keys() == set(['Berger_POPC', 'Berger_PLA', 'Berger_POP', 'CHARMM_POPC'])

        # Test if each dic contains the right number of atoms
        for name in lipids_tops.keys():
            assert len(lipids_tops[name]) == reference[name]


    def test_check_read_lipids_topH_failure(self):
        """Test check_read_lipids_topH() with a file in a wrong format."""
        bad_file = "Berger_wrongformat.json"
        with pytest.raises(UserWarning) as err:
            lipids.read_lipids_topH([path_data / bad_file])
        assert f"{bad_file} is in a bad format." in str(err.value)
