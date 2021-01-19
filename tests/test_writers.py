"""
Unit tests for buildH_calcOP.

Test functions from module writers.
"""

import pathlib
import filecmp
import pytest

import pandas as pd
import MDAnalysis as mda

from buildh import lipids
from buildh import init_dics
from buildh import core
from buildh import writers


DIR_DATA = "test_data"
PATH_ROOT_DATA = pathlib.Path(__file__).parent / DIR_DATA

# Ignore some MDAnalysis warnings for this test file
pytestmark = pytest.mark.filterwarnings('ignore::UserWarning')


def test_pandas2pdb():
    """Test for pandasdf2pdb()."""
    # Create a dummy dataframe
    rows = [[1, "C1",   "POPC", 1, 34.42, 46.94, 26.31],
            [2, "H211", "POPC", 1,  1.00,  2.00,  3.00]
    ]
    df = pd.DataFrame(rows, columns=["atnum", "atname", "resname", "resnum",
                                     "x", "y", "z"])

    ref_pdb_lines = ('ATOM      1  C1  POPC    1      34.420  46.940  26.310  1.00  0.00             C\n'
                     'ATOM      2 H211 POPC    1       1.000   2.000   3.000  1.00  0.00             H\n')

    assert writers.pandasdf2pdb(df) == ref_pdb_lines


class TestWriters():
    """Test class for the writing functions."""

    # path for the Berger POPC files
    PATH_DATA = PATH_ROOT_DATA / "Berger_POPC"

    # Method called once per class.
    def setup_class(self):
        """Initialize all data."""
        lipids_tops = lipids.read_lipids_topH([lipids.PATH_JSON/"Berger_POPC.json"])

        # Input parameters
        self.pdb = self.PATH_DATA / "10POPC.pdb"
        self.defop = self.PATH_DATA / "OP_def_BergerPOPC.def"
        self.dic_lipid = lipids_tops["Berger_POPC"]

        self.begin = 0
        self.end = 1

        # attributes
        self.universe_woH = mda.Universe(str(self.pdb))
        self.dic_atname2genericname = init_dics.make_dic_atname2genericname(self.defop)
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid['resname'])
        self.dic_Cname2Hnames = init_dics.make_dic_Cname2Hnames(self.dic_OP)

        # Compute the order parameter
        core.fast_build_all_Hs_calc_OP(self.universe_woH,self.begin, self.end,
                                       self.dic_OP, self.dic_lipid, self.dic_Cname2Hnames)

    def test_write_OP(self, tmpdir):
        """
        Test for write_OP().

        Parameters
        ----------
        tmpdir: function
            pytest callback which return a unique directory.
        """
        #Write results of self.dic_OP
        test_file = tmpdir / "test.out"
        writers.write_OP(test_file, self.dic_atname2genericname,
                                self.dic_OP, self.dic_lipid['resname'])

        ref_file = self.PATH_DATA / "ref_10POPC.out"
        assert filecmp.cmp(test_file, ref_file)

    def test_write_OP_alternate(self, tmpdir):
        """
        Test for write_OP_alternate().

        Parameters
        ----------
        tmpdir: function
            pytest callback which return a unique directory.
        """
        #Write results of self.dic_OP
        test_file = tmpdir / "test.out"
        writers.write_OP_alternate(test_file, self.universe_woH,
                                   self.dic_OP, self.dic_lipid['resname'])

        ref_file = self.PATH_DATA / "ref_10POPC.alternate.out"
        assert filecmp.cmp(test_file, ref_file)
