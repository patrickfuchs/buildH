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
import pandas as pd

from buildh import dic_lipids
from buildh import init_dics
from buildh import core
from buildh import writers

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
        self.universe_woH = mda.Universe(str(self.pdb))
        self.dic_atname2genericname = init_dics.make_dic_atname2genericname(self.defop)
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid)
        self.dic_Cname2Hnames = init_dics.make_dic_Cname2Hnames(self.dic_OP)


    # Method called before each test method.
    def setup_method(self):
        """
        self.dic_op needs to be reinitialized since it is modified by some of the functions tested.
        """
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid)

    ########################################
    # Tests methods for the fast algorithm #
    ########################################

    @pytest.mark.parametrize('atom_index, helpers_indexes', [
         # atom C1
         (0, (3,4)),
         # atom C13
         (12, (11, 31, 13)),
         # atom C24
         (23, (22, 24)),
         # atom C50
         (49, (48, 47)),
    ])
    def test_get_indexes(self, atom_index, helpers_indexes):
        """
        Test for get_indexes()
        """
        # Use the value in parametrize to get the correct atom to test.
        atom = self.universe_woH.atoms[atom_index]

        assert core.get_indexes(atom, self.dic_lipid) == helpers_indexes


    def test_make_dic_lipids_with_indexes(self):
        """
        Test for make_dic_lipids_with_indexes
        """

        dic_lipids_with_indexes = core.make_dic_lipids_with_indexes(self.universe_woH,
                                                                     self.dic_lipid,
                                                                     self.dic_OP)

        assert dic_lipids_with_indexes['C1']  == ('CH3', 'N4', 'C5', 0, 3, 4)
        assert dic_lipids_with_indexes['C5']  == ('CH2', 'N4', 'C6', 4, 3, 5)
        assert dic_lipids_with_indexes['C13'] == ('CH', 'C12', 'C32', 'O14', 12, 11, 31, 13)
        assert dic_lipids_with_indexes['C24'] == ('CHdoublebond', 'C23', 'C25', 23, 22, 24)
        assert dic_lipids_with_indexes['C50'] == ('CH3', 'C49', 'C48', 49, 48, 47)


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
        The results should be indentical to the test_reconstruct_Hs() test.

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

    ########################################
    # Tests methods for the slow algorithm #
    # where hydrogens are saved            #
    ########################################

    # References values are the same as in test_fast_buildHs_on_1C()
    @pytest.mark.parametrize('index, Hs_coords', [
            # atom C1 type CH3
            (0,  (np.array([35.06161421, 47.69320272, 26.76728762]),
                  np.array([33.93128850, 47.36328732, 25.43245201]),
                  np.array([35.02249090, 46.08195972, 26.01188584]))),
            # atom C5 type CH2
            (4,  (np.array([31.62418725, 45.50474260, 27.61469385]),
                  np.array([32.93474471, 44.52992542, 26.90727792]))),
            # atom C13 type CH
            (12, (np.array([26.61868981, 44.14296091, 25.55244500]),)),
            # atom C24 type CHdoublebond
            (23, (np.array([20.46454439, 38.99866187, 31.19224364]),)),
    ])
    def test_buildHs_on_1C(self, index, Hs_coords):
        """
        Test for buildHs_on_1C()
        Generate 4 atoms to be tested, each with a different type
        """
        atom = self.universe_woH.atoms[index]
        test_Hs_coords = core.buildHs_on_1C(atom, self.dic_lipid)

        assert_almost_equal(test_Hs_coords, Hs_coords)


    def test_reconstruct_Hs_first_frame(self, tmpdir):
        """
        Test for build_all_Hs_calc_OP() in the first mode
        """
        new_df_atoms = core.build_all_Hs_calc_OP(self.universe_woH, self.dic_lipid,
                                                 self.dic_Cname2Hnames, return_coors=True)

        assert new_df_atoms.shape == (1340, 7)

        # Test random atoms
        col_series=["atnum", "atname","resname", "resnum","x", "y", "z"]
        ref_atom = pd.Series([1, "C1", "POPC", 2, 34.41999816, 46.93999862, 26.30999946], index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[0], ref_atom, check_names=False)

        ref_atom = pd.Series([2, "H11", "POPC", 2, 35.06161421, 47.69320272, 26.76728762], index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[1], ref_atom, check_names=False)

        ref_atom = pd.Series([259, "H501", "POPC", 3, 74.16520229, 36.21701104, 35.03256486], index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[258], ref_atom, check_names=False)

        ref_atom = pd.Series([1339, "HA22", "POPC", 11, 27.72946942, 16.74704078, 40.53260384], index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[1338], ref_atom, check_names=False)


    def test_reconstruct_Hs(self, tmpdir):
        """
        Test for build_all_Hs_calc_OP() in the second mode
        The results should be indentical to the test_fast_build_all_Hs_calc_OP() test.
        """
        pdb_wH = path_data / "10POPC_wH.pdb"
        universe_wH = mda.Universe(str(pdb_wH))
        core.build_all_Hs_calc_OP(self.universe_woH, self.dic_lipid, self.dic_Cname2Hnames,
                             universe_wH=universe_wH, dic_OP=self.dic_OP,
                             dic_corresp_numres_index_dic_OP=self.dic_corresp_numres_index_dic_OP)


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



# Ignore some MDAnalysis warnings
@pytest.mark.filterwarnings('ignore::UserWarning')
class TestXTCPOPC:
    """
    Test class for a trajectory of POPC lipids.
    """

    # Method called once per class.
    def setup_class(self):
        """
        Initialize all data.
        """

        # Input parameters
        self.pdb = path_data / "2POPC.pdb"
        self.xtc = path_data / "2POPC.xtc"
        self.defop = path_data / "OP_def_BergerPOPC.def"
        self.dic_lipid = getattr(dic_lipids, "Berger_POPC")
        self.begin = 0
        self.end = 11

        # attributes
        self.universe_woH = mda.Universe(str(self.pdb), str(self.xtc))
        self.dic_atname2genericname = init_dics.make_dic_atname2genericname(self.defop)
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid)
        self.dic_Cname2Hnames = init_dics.make_dic_Cname2Hnames(self.dic_OP)

    # Method called before each test method.
    def setup_method(self):
        """
        self.dic_op needs to be reinitialized since it is modified by some of the functions tested.
        """
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid)

    def test_fast_calcOP(self, tmpdir):
        """
        Test fast_build_all_Hs_calc_OP() on a trajectory
        The results should be indentical to the test_gen_XTC_calcOP() test.
        """

        core.fast_build_all_Hs_calc_OP(self.universe_woH,self.begin, self.end,
                                       self.dic_OP, self.dic_lipid, self.dic_Cname2Hnames)

        test_file_jmelcr = tmpdir / "test_2POPC_traj.jmelcr.out"
        test_file_apineiro = tmpdir / "test_2POPC_traj.apineiro.out"
        writers.write_OP_jmelcr(test_file_jmelcr, self.dic_atname2genericname,
                                self.dic_OP, self.dic_lipid)
        writers.write_OP_apineiro(test_file_apineiro, self.universe_woH,
                                  self.dic_OP, self.dic_lipid)

        ref_file_jmelcr = path_data / "ref_2POPC_traj.jmelcr.out"
        ref_file_apineiro = path_data / "ref_2POPC_traj.apineiro.out"
        assert filecmp.cmp(test_file_jmelcr, ref_file_jmelcr)
        assert filecmp.cmp(test_file_apineiro, ref_file_apineiro)

    def test_gen_XTC_calcOP(self, tmpdir):
        """
        Test for gen_XTC_calcOP()
        The results should be indentical to the test_fast_calcOP() test.
        """
        core.gen_XTC_calcOP("test", self.universe_woH, self.dic_OP, self.dic_lipid,
                            self.dic_Cname2Hnames, self.dic_corresp_numres_index_dic_OP,
                            self.begin, self.end)

        #Write results
        test_file_jmelcr = tmpdir / "test_2POPC_traj.jmelcr.out"
        test_file_apineiro = tmpdir / "test_2POPC_traj.apineiro.out"
        writers.write_OP_jmelcr(test_file_jmelcr, self.dic_atname2genericname,
                                self.dic_OP, self.dic_lipid)
        writers.write_OP_apineiro(test_file_apineiro, self.universe_woH,
                                  self.dic_OP, self.dic_lipid)

        ref_file_jmelcr = path_data / "ref_2POPC_traj.jmelcr.out"
        ref_file_apineiro = path_data / "ref_2POPC_traj.apineiro.out"
        assert filecmp.cmp(test_file_jmelcr, ref_file_jmelcr)
        assert filecmp.cmp(test_file_apineiro, ref_file_apineiro)
