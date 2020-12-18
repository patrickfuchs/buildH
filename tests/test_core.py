#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""
Unit tests for buildH_calcOP

Test functions from module core
"""

import pathlib
import filecmp
import pytest

import numpy as np
from numpy.testing import assert_almost_equal
import MDAnalysis as mda

import dic_lipids
import init_dics
import core
import writers

dir_data = "test_data"
path_data = pathlib.Path(__file__).parent / dir_data

@pytest.fixture(scope='class')
def tmp_dir(tmp_path_factory):
    """
    Create a temp directory for the result files
    """
    return tmp_path_factory.mktemp("data")

# Ignore some MDAnalysis warnings
@pytest.mark.filterwarnings('ignore::UserWarning')
class TestPDBPOPC:
    """
    Test class for a pdb file (no traj) of POPC lipids.
    """

    # Method called once per class.
    def setup_class(self):
        """
        Initialize all data.
        """

        # Input parameters
        self.pdb = path_data / "10POPC.pdb"
        self.defop = path_data / "OP_def_BergerPOPC.def"
        self.dic_lipid = getattr(dic_lipids, "Berger_POPC")
        self.begin = 0
        self.end = 1

        # attributes
        self.universe_woH = mda.Universe(self.pdb)
        self.dic_atname2genericname = init_dics.make_dic_atname2genericname(self.defop)
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid)
        self.dic_Cname2Hnames = init_dics.make_dic_Cname2Hnames(self.dic_OP)


    # Method called before each test method.
    def setup_method(self):
        """
        self.dic_op needs to be reinitialized since it is modified by some of the functions tested.

        TODO: For now, only one function modify it. Maybe move this in setup_class() if it remains like this.
        """
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid)


    @pytest.mark.parametrize('Cname, index, Hs_coords', [
            # type CH3
            ("C1", 0, (np.array([35.06161421, 47.69320272, 26.76728762]),
                       np.array([33.93128850, 47.36328732, 25.43245201]),
                       np.array([35.02249090, 46.08195972, 26.01188584]))),
            # type CH2
            ("C5", 0, (np.array([31.62418725, 45.50474260, 27.61469385]),
                       np.array([32.93474471, 44.52992542, 26.90727792]))),
            # type CH
            ("C13", 0, (np.array([26.61868981, 44.14296091, 25.55244500]),)),
            # type CHdoublebond
            ("C24", 0, (np.array([20.46454439, 38.99866187, 31.19224364]),)),
    ])
    def test_fast_buildHs_on_1C(self, Cname, index, Hs_coords):
        """
            Test for fast_buildHs_on_1C()
            Generate 4 atoms to be tested, each with a different type
        """

        dic_lipids_with_indexes = core.make_dic_lipids_with_indexes(self.universe_woH,
                                                                     self.dic_lipid,
                                                                     self.dic_OP)
        ts = self.universe_woH.trajectory[0]
        test_Hs_coords = core.fast_buildHs_on_1C(dic_lipids_with_indexes, ts, Cname, index)

        assert_almost_equal(test_Hs_coords, Hs_coords)


    def test_fast_build_all_Hs_calc_OP(self, tmpdir):
        """
            Test for fast_build_all_Hs_calc_OP()

            TODO: For now, we test the result file so we test different functions in one test :
            fast_build_all_Hs_calc_OP and write functions
            It should be splitted.
        """

        core.fast_build_all_Hs_calc_OP(self.universe_woH,self.begin, self.end,
                                       self.dic_OP, self.dic_lipid, self.dic_Cname2Hnames)

        #Write results
        test_file_jmelcr = tmpdir / "test_10POPC.jmelcr.out"
        test_file_apineiro = tmpdir / "test_10POPC.apineiro.out"
        writers.write_OP_jmelcr(test_file_jmelcr, self.dic_atname2genericname,
                                self.dic_OP, self.dic_lipid)
        writers.write_OP_apineiro(test_file_apineiro, self.universe_woH,
                                  self.dic_OP, self.dic_lipid)

        ref_file_jmelcr = path_data / "ref_10POPC.jmelcr.out"
        ref_file_apineiro = path_data / "ref_10POPC.apineiro.out"
        assert filecmp.cmp(test_file_jmelcr, ref_file_jmelcr)
        assert filecmp.cmp(test_file_apineiro, ref_file_apineiro)



