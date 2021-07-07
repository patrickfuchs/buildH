"""
Unit tests for buildH_calcOP.

Test functions from module core
"""

import os
import pathlib
import pytest

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
import MDAnalysis as mda
import pandas as pd

from buildh import lipids
from buildh import init_dics
from buildh import core

DIR_DATA = "test_data"
PATH_ROOT_DATA = pathlib.Path(__file__).parent / DIR_DATA


# Ignore some MDAnalysis warnings
@pytest.mark.filterwarnings('ignore::UserWarning')
class TestPDBPOPC:
    """Test class for a single pdb file of POPC lipids."""

    # path for the Berger POPC files
    PATH_DATA = PATH_ROOT_DATA / "Berger_POPC"

    # Subset of reference data for dic_OP result
    # This used for test_fast_build_all_Hs_calc_OP() and test_reconstruct_Hs()
    ref_OP = {
        # values = 10 lipids x 1 frame
        ('C1', 'H11'):   [[-0.2359911], [ 0.5787392], [ 0.0224484], [ 0.32463246], [-0.39586149],
                          [ 0.87975540], [ 0.41808751], [-0.28873301], [-0.47090728], [-0.11264333]],

        ('C32', 'H321'): [[-0.41362472], [-0.49548676], [ 0.14942140], [ 0.94630521], [ 0.31010481],
                          [ 0.89795911], [-0.49998903], [-0.49960781], [ 0.79312354], [-0.05705253]],

        ('C19', 'H192'): [[-0.00139755], [-0.26007410], [-0.42209355], [-0.27136414], [-0.35669670],
                          [-0.46162091], [-0.25302619], [ 0.65296792], [ 0.14922587], [-0.41107551]],

        ('CA1', 'HA11'): [[-0.48614851], [-0.43749245], [-0.40528734], [-0.35136805], [-0.45665044],
                          [-0.21534620], [ 0.08544723], [ 0.44870036], [-0.10571545], [-0.46312000]]
    }

    # Method called once per class.
    def setup_class(self):
        """Initialize attributes."""
        lipids_tops = lipids.read_lipids_topH([lipids.PATH_JSON/"Berger_POPC.json"])

        # Input parameters
        self.pdb = self.PATH_DATA / "10POPC.pdb"
        self.defop = self.PATH_DATA / "OP_def_BergerPOPC.def"
        self.dic_lipid = lipids_tops["Berger_POPC"]
        self.begin = 0
        self.end = 1

        # intern attributes
        self.universe_woH = mda.Universe(str(self.pdb))
        self.dic_atname2genericname = init_dics.make_dic_atname2genericname(self.defop)
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid['resname'])
        self.dic_Cname2Hnames = init_dics.make_dic_Cname2Hnames(self.dic_OP)


    # Method called before each test method.
    def setup_method(self):
        """Reset some attributes.

        self.dic_op needs to be reinitialized
        since it is modified by some of the functions tested.
        """
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid['resname'])


    @pytest.mark.parametrize('atom_index, Hs_coords', [
            # atom C1 type CH3
            (0,  (np.array([35.06161421, 47.69320272, 26.76728762], dtype=np.float32),
                  np.array([33.93128850, 47.36328732, 25.432453], dtype=np.float32),
                  np.array([35.02249090, 46.08195972, 26.011885], dtype=np.float32))),
            # atom C5 type CH2
            (4,  (np.array([31.62418725, 45.50474260, 27.61469385], dtype=np.float32),
                  np.array([32.93474471, 44.52992542, 26.90727792], dtype=np.float32))),
            # atom C13 type CH
            (12, (np.array([26.6186898,  44.1429600, 25.552444], dtype=np.float32),)),
            # atom C24 type CHdoublebond
            (23, (np.array([20.46454439, 38.99866187, 31.19224364], dtype=np.float32),)),
    ])
    def test_buildHs_on_1C(self, atom_index, Hs_coords):
        """Test for buildHs_on_1C().

        Generate 4 atoms to be tested, each with a different type.

        Parameters
        ----------
        atom_index : int
            index of the heavy atom where hydrogens will be reconstructed.
        Hs_coords : numpy 1D-array
            reference coordinates of the hydrogens rebuilt.
        """
        atom = self.universe_woH.atoms[atom_index]

        sel = atom.residue.atoms.select_atoms
        if len(self.dic_lipid[atom.name]) == 3:
            typeofH2build, helper1_name, helper2_name = self.dic_lipid[atom.name]
            helper3_coor = None
        else:
            typeofH2build, helper1_name, helper2_name, helper3_name = self.dic_lipid[atom.name]
            helper3_coor = sel("name {0}".format(helper3_name))[0].position

        helper1_coor = sel("name {0}".format(helper1_name))[0].position
        helper2_coor = sel("name {0}".format(helper2_name))[0].position


        test_Hs_coords = core.buildHs_on_1C(atom.position, typeofH2build,
                                            helper1_coor, helper2_coor, helper3_coor)

        assert_almost_equal(test_Hs_coords, Hs_coords)


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
        """Test for get_indexes().

        Parameters
        ----------
        atom_index : int
            index of a heavy atom
        helpers_indexes : tuple
            references indexes of the helpers
        """
        # Use the value in parametrize to get the correct atom to test.
        atom = self.universe_woH.atoms[atom_index]

        assert core.get_indexes(atom, self.dic_lipid) == helpers_indexes


    def test_make_dic_lipids_with_indexes(self):
        """Test for make_dic_lipids_with_indexes()."""
        dic_lipids_with_indexes = core.make_dic_lipids_with_indexes(self.universe_woH,
                                                                     self.dic_lipid,
                                                                     self.dic_OP)

        assert dic_lipids_with_indexes['C1']  == ['CH3', 'N4', 'C5', 0, 3, 4]
        assert dic_lipids_with_indexes['C5']  == ['CH2', 'N4', 'C6', 4, 3, 5]
        assert dic_lipids_with_indexes['C13'] == ['CH', 'C12', 'C32', 'O14', 12, 11, 31, 13]
        assert dic_lipids_with_indexes['C24'] == ['CHdoublebond', 'C23', 'C25', 23, 22, 24]
        assert dic_lipids_with_indexes['C50'] == ['CH3', 'C49', 'C48', 49, 48, 47]


    def test_fast_build_all_Hs_calc_OP(self):
        """Test for fast_build_all_Hs_calc_OP().

        The results should be indentical to the test_reconstruct_Hs() test.
        """
        core.fast_build_all_Hs_calc_OP(self.universe_woH,self.begin, self.end,
                                       self.dic_OP, self.dic_lipid, self.dic_Cname2Hnames)

        # Check statistics
        assert_almost_equal(np.mean(self.dic_OP[('C40', 'H401')]), -0.28794656)
        assert_almost_equal(np.mean(self.dic_OP[('C17', 'H171')]), -0.18843357)

        # Check few particular cases
        # Use of a loop to check key and value separately.
        # The values need to be assert with float desired precision.
        for (key), value in self.ref_OP.items():
            assert key in self.dic_OP.keys()
            assert_almost_equal(self.dic_OP[key], value)



    ########################################
    # Tests methods for the slow algorithm #
    # where hydrogens are saved            #
    ########################################


    def test_reconstruct_Hs_first_frame(self):
        """Test for build_system_hydrogens()."""
        dic_lipids_with_indexes = core.make_dic_lipids_with_indexes(self.universe_woH,
                                                                     self.dic_lipid,
                                                                     self.dic_OP)

        new_df_atoms = core.build_system_hydrogens(self.universe_woH, self.dic_lipid,
                                                   self.dic_Cname2Hnames, dic_lipids_with_indexes)

        assert new_df_atoms.shape == (1340, 7)

        # Test random atoms
        col_series=["atnum", "atname","resname", "resnum","x", "y", "z"]
        ref_atom = pd.Series([1, "C1", "POPC", 2, 34.41999816, 46.93999862, 26.30999946],
                             index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[0], ref_atom, check_names=False)

        ref_atom = pd.Series([2, "H11", "POPC", 2, 35.06161421, 47.69320272, 26.76728762],
                             index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[1], ref_atom, check_names=False)

        ref_atom = pd.Series([259, "H501", "POPC", 3, 74.16520229, 36.21701104, 35.03256486],
                             index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[258], ref_atom, check_names=False)

        ref_atom = pd.Series([1339, "HA22", "POPC", 11, 27.72946942, 16.74704078, 40.53260384],
                             index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[1338], ref_atom, check_names=False)


    def test_reconstruct_Hs(self):
        """Test for build_all_Hs_calc_OP().

        The results should be indentical to the test_fast_build_all_Hs_calc_OP() test.
        """
        dic_lipids_with_indexes = core.make_dic_lipids_with_indexes(self.universe_woH,
                                                                    self.dic_lipid,
                                                                    self.dic_OP)

        pdb_wH = self.PATH_DATA / "10POPC_wH.pdb"
        universe_wH = mda.Universe(str(pdb_wH))
        ts = self.universe_woH.trajectory[0]
        core.build_all_Hs_calc_OP(self.universe_woH, ts, self.dic_lipid, self.dic_Cname2Hnames,
                                  universe_wH, self.dic_OP, self.dic_corresp_numres_index_dic_OP,
                                  dic_lipids_with_indexes)

        # Check statistics
        assert_almost_equal(np.mean(self.dic_OP[('C21', 'H212')]), -0.22229490)
        assert_almost_equal(np.mean(self.dic_OP[('CA2', 'HA23')]),  0.22305829)

        # Check few particular cases
        # Use of a loop to check key and value separately.
        # The values need to be assert with float desired precision.
        for (key), value in self.ref_OP.items():
            assert key in self.dic_OP.keys()
            assert_almost_equal(value, self.dic_OP[key])



# Ignore some MDAnalysis warnings
@pytest.mark.filterwarnings('ignore::UserWarning')
class TestXTCPOPC:
    """Test class for a trajectory (in xtc format) of POPC lipids."""

    # path for the Berger POPC files
    PATH_DATA = PATH_ROOT_DATA / "Berger_POPC"

    # Subset of reference data for dic_OP result
    # This used for test_fast_calcOP() and test_gen_coordinates_calcOP()
    ref_OP = {
        # values = 2 lipids x 11 frames
        ('C2', 'H23'):   [ [-0.4729392, -0.4853177, -0.3329995, -0.0827970, -0.4993968,
                            -0.4990348, -0.3598688,  0.7689081,  0.0815748,  0.1942028,
                             0.2798647],
                           [ 0.1614778, -0.3463071,  0.4329743,  0.2873318, -0.2959748,
                            -0.3995778, -0.4891822, -0.0830738,  0.5635600,  0.7569330,
                             0.2006300]
                         ],
        ('C27', 'H272'): [ [ 0.2147729, -0.4988779, -0.4161010, -0.1932931,  0.3062109,
                             0.5009914, -0.2168313,  0.5028035,  0.0603169, -0.0039737,
                            -0.4971311],
                           [ 0.0128815, -0.4654747, -0.3786469, -0.2972043, -0.4199277,
                            -0.0001302, -0.3428786, -0.4642762, -0.2647048,  0.8727489,
                            -0.2021679]
                         ]
    }

    # Method called once per class.
    def setup_class(self):
        """Initialize attributes."""
        lipids_tops = lipids.read_lipids_topH([lipids.PATH_JSON/"Berger_POPC.json"])

        # Input parameters
        self.pdb = self.PATH_DATA / "2POPC.pdb"
        self.xtc = self.PATH_DATA / "2POPC.xtc"
        self.defop = self.PATH_DATA / "OP_def_BergerPOPC.def"
        self.dic_lipid = lipids_tops["Berger_POPC"]
        self.begin = 0
        self.end = 11

        # attributes
        self.universe_woH = mda.Universe(str(self.pdb), str(self.xtc))
        self.dic_atname2genericname = init_dics.make_dic_atname2genericname(self.defop)
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid['resname'])
        self.dic_Cname2Hnames = init_dics.make_dic_Cname2Hnames(self.dic_OP)

    # Method called before each test method.
    def setup_method(self):
        """Reset some attributes.

        self.dic_op needs to be reinitialized
        since it is modified by some of the functions tested.
        """
        self.dic_OP, self.dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(self.universe_woH,
                                                                                  self.dic_atname2genericname,
                                                                                  self.dic_lipid['resname'])

    def test_fast_calcOP(self):
        """Test fast_build_all_Hs_calc_OP() on a trajectory.

        The results should be indentical to the test_gen_coordinates_calcOP() test.
        """
        core.fast_build_all_Hs_calc_OP(self.universe_woH,self.begin, self.end,
                                       self.dic_OP, self.dic_lipid, self.dic_Cname2Hnames)

        # Check statistics
        assert_almost_equal(np.mean(self.dic_OP[('C32', 'H321')]),   0.1530018)
        assert_almost_equal(np.mean(self.dic_OP[('C50', 'H503')]), -0.08801112)
        assert_almost_equal(np.mean(self.dic_OP[('C1', 'H11')]),  0.26908040)
        assert_almost_equal(np.mean(self.dic_OP[('C5', 'H52')]), -0.20147210)

        # Check few particular cases
        # Use of a loop to check key and value separately.
        # The values need to be assert with float desired precision.
        for (key), value in self.ref_OP.items():
            assert key in self.dic_OP.keys()
            assert_almost_equal(self.dic_OP[key], value)

    def test_gen_coordinates_calcOP(self):
        """Test for gen_coordinates_calcOP().

        The results should be indentical to the test_fast_calcOP() test.
        """
        core.gen_coordinates_calcOP("test", self.universe_woH, self.dic_OP, self.dic_lipid,
                                    self.dic_Cname2Hnames, self.dic_corresp_numres_index_dic_OP,
                                    self.begin, self.end, True)

        # Check  statistics
        assert_almost_equal(np.mean(self.dic_OP[('C32', 'H321')]),  0.15300180)
        assert_almost_equal(np.mean(self.dic_OP[('C50', 'H503')]), -0.08801112)
        assert_almost_equal(np.mean(self.dic_OP[('C1', 'H11')]),  0.26908040)
        assert_almost_equal(np.mean(self.dic_OP[('C5', 'H52')]), -0.20147210)

        # Check few particular cases
        # Use of a loop to check key and value separately.
        # The values need to be assert with float desired precision.
        for (key), value in self.ref_OP.items():
            assert key in self.dic_OP.keys()
            assert_almost_equal(self.dic_OP[key], value)


    def test_output_traj(self, tmp_path):
        """Test the validity of the ouput trajectory files when requested

        Parameters
        ----------
        tmp_path: function
            pytest callback which return a unique directory.
        """

        os.chdir(tmp_path)
        core.gen_coordinates_calcOP("test", self.universe_woH, self.dic_OP, self.dic_lipid,
                                    self.dic_Cname2Hnames, self.dic_corresp_numres_index_dic_OP,
                                    self.begin, self.end, True)

        #Tests macro values of the output files
        u = mda.Universe("test.pdb", "test.xtc")
        assert_equal(u.trajectory.n_frames, 11)
        assert_equal(u.atoms.n_atoms, 268)
        assert_almost_equal(u.trajectory[4].dimensions,
                            np.array([66.251070, 66.25107, 88.025925, 90.0, 90.0, 90.0], dtype=np.float32))

