"""
Unit tests for buildH_calcOP.

Test functions from module core
"""

import pathlib
import pytest

import numpy as np
from numpy.testing import assert_almost_equal
import MDAnalysis as mda
import pandas as pd

from buildh import lipids
from buildh import init_dics
from buildh import core
from buildh import writers

dir_data = "test_data"
path_data = pathlib.Path(__file__).parent / dir_data


# Ignore some MDAnalysis warnings
@pytest.mark.filterwarnings('ignore::UserWarning')
class TestPDBPOPC:
    """Test class for a single pdb file of POPC lipids."""

    # Subset of reference data for dic_OP result
    # This used for test_fast_build_all_Hs_calc_OP() and test_reconstruct_Hs()
    ref_OP = {
            # values = 10 lipids x 1 frame
            ('C1', 'H11'):   [[-0.23599134], [ 0.57873978], [ 0.02244849], [ 0.32463168], [-0.39586207],
                              [ 0.87975500], [ 0.41808654], [-0.28873204], [-0.47090717], [-0.11264325]],

            ('C32', 'H321'): [[-0.41362438], [-0.49548668], [ 0.14942125], [ 0.94630508], [ 0.31010313],
                              [ 0.89795921], [-0.49998903], [-0.49960781], [ 0.79312234], [-0.05705310]],

            ('C19', 'H192'): [[-0.00139740], [-0.26007547], [-0.42209310], [-0.27136435], [-0.35669634],
                              [-0.46162108], [-0.25302670], [ 0.65296859], [ 0.14922439], [-0.41107602]],

            ('CA1', 'HA11'): [[-0.48614830], [-0.43749174], [-0.40528699], [-0.35136869], [-0.45664981],
                              [-0.21534578], [ 0.08544719], [ 0.44870156], [-0.10571712], [-0.46312028]]
    }

    # Method called once per class.
    def setup_class(self):
        """Initialize attributes."""
        lipids_tops = lipids.read_lipids_topH([lipids.PATH_JSON/"Berger_POPC.json"])

        # Input parameters
        self.pdb = path_data / "10POPC.pdb"
        self.defop = path_data / "OP_def_BergerPOPC.def"
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


        test_Hs_coords = core.buildHs_on_1C(atom.position, typeofH2build, helper1_coor, helper2_coor, helper3_coor)

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
        # Use of a loop to check key and value separately. The values need to be assert with float desired precision.
        for (key), value in self.ref_OP.items():
            assert key in self.dic_OP.keys()
            assert_almost_equal(value, self.dic_OP[key])


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
        ref_atom = pd.Series([1, "C1", "POPC", 2, 34.41999816, 46.93999862, 26.30999946], index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[0], ref_atom, check_names=False)

        ref_atom = pd.Series([2, "H11", "POPC", 2, 35.06161421, 47.69320272, 26.76728762], index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[1], ref_atom, check_names=False)

        ref_atom = pd.Series([259, "H501", "POPC", 3, 74.16520229, 36.21701104, 35.03256486], index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[258], ref_atom, check_names=False)

        ref_atom = pd.Series([1339, "HA22", "POPC", 11, 27.72946942, 16.74704078, 40.53260384], index=col_series)
        pd.testing.assert_series_equal(new_df_atoms.loc[1338], ref_atom, check_names=False)


    def test_reconstruct_Hs(self):
        """Test for build_all_Hs_calc_OP().

        The results should be indentical to the test_fast_build_all_Hs_calc_OP() test.
        """
        dic_lipids_with_indexes = core.make_dic_lipids_with_indexes(self.universe_woH,
                                                                    self.dic_lipid,
                                                                    self.dic_OP)

        pdb_wH = path_data / "10POPC_wH.pdb"
        universe_wH = mda.Universe(str(pdb_wH))
        ts = self.universe_woH.trajectory[0]
        core.build_all_Hs_calc_OP(self.universe_woH, ts, self.dic_lipid, self.dic_Cname2Hnames,
                                  universe_wH, self.dic_OP, self.dic_corresp_numres_index_dic_OP,
                                  dic_lipids_with_indexes)

        # Check statistics
        assert_almost_equal(np.mean(self.dic_OP[('C21', 'H212')]), -0.22229490)
        assert_almost_equal(np.mean(self.dic_OP[('CA2', 'HA23')]),  0.22305860)

        # Check few particular cases
        # Use of a loop to check key and value separately. The values need to be assert with float desired precision.
        for (key), value in self.ref_OP.items():
            assert key in self.dic_OP.keys()
            assert_almost_equal(value, self.dic_OP[key])



# Ignore some MDAnalysis warnings
@pytest.mark.filterwarnings('ignore::UserWarning')
class TestXTCPOPC:
    """Test class for a trajectory (in xtc format) of POPC lipids."""

    # Subset of reference data for dic_OP result
    # This used for test_fast_calcOP() and test_gen_coordinates_calcOP()
    ref_OP = {
        # values = 2 lipids x 11 frames
        ('C2', 'H23'):   [ [-0.47293904, -0.48531776, -0.33300023, -0.08279667, -0.49939686,
                            -0.49903488, -0.35986740,  0.76890823,  0.08157466,  0.19420324, 0.27986448],
                           [ 0.16147857, -0.34630697,  0.43297339,  0.28733088, -0.29597471,
                            -0.39957748, -0.48918187, -0.08307259,  0.56356016,  0.75693095,  0.20062960]
                         ],
        ('C27', 'H272'): [ [ 0.21477296, -0.49887799, -0.41609987, -0.19329433,  0.30621276,
                             0.50099123, -0.21682994,  0.50280378,  0.06031779, -0.00397443, -0.49713104],
                           [ 0.01288211, -0.46547434, -0.37864825, -0.29720486, -0.41992799,
                            -0.00013077, -0.342879725, -0.464275827, -0.264702934, 0.87274949,-0.20216758]
                         ]
    }

    # Method called once per class.
    def setup_class(self):
        """Initialize attributes."""
        lipids_tops = lipids.read_lipids_topH([lipids.PATH_JSON/"Berger_POPC.json"])

        # Input parameters
        self.pdb = path_data / "2POPC.pdb"
        self.xtc = path_data / "2POPC.xtc"
        self.defop = path_data / "OP_def_BergerPOPC.def"
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

    def test_fast_calcOP(self, tmp_path):
        """Test fast_build_all_Hs_calc_OP() on a trajectory.

        The results should be indentical to the test_gen_coordinates_calcOP() test.

        Parameters
        ----------
        tmp_path : pathlib.Path (Pytest fixture)
            path to a unique temporary directory.
        """
        core.fast_build_all_Hs_calc_OP(self.universe_woH,self.begin, self.end,
                                       self.dic_OP, self.dic_lipid, self.dic_Cname2Hnames)

        # Check statistics
        assert_almost_equal(np.mean(self.dic_OP[('C1', 'H11')]),  0.26908040)
        assert_almost_equal(np.mean(self.dic_OP[('C5', 'H52')]), -0.20147210)

        # Check few particular cases
        # Use of a loop to check key and value separately. The values need to be assert with float desired precision.
        for (key), value in self.ref_OP.items():
            assert key in self.dic_OP.keys()
            assert_almost_equal(value, self.dic_OP[key])

    def test_gen_coordinates_calcOP(self, tmp_path):
        """Test for gen_coordinates_calcOP().

        The results should be indentical to the test_fast_calcOP() test.

        Parameters
        ----------
        tmp_path : pathlib.Path (Pytest fixture)
            path to a unique temporary directory.
        """
        core.gen_coordinates_calcOP("test", self.universe_woH, self.dic_OP, self.dic_lipid,
                                    self.dic_Cname2Hnames, self.dic_corresp_numres_index_dic_OP,
                                    self.begin, self.end, True)

        # Check  statistics
        assert_almost_equal(np.mean(self.dic_OP[('C32', 'H321')]),  0.15300163)
        assert_almost_equal(np.mean(self.dic_OP[('C50', 'H503')]), -0.08801085)

        # Check few particular cases
        # Use of a loop to check key and value separately. The values need to be assert with float desired precision.
        for (key), value in self.ref_OP.items():
            assert key in self.dic_OP.keys()
            assert_almost_equal(value, self.dic_OP[key])


    def test_check_def_file(self):
        """Test for is_allHs_present()."""
        assert core.is_allHs_present(self.defop, self.dic_lipid, self.dic_Cname2Hnames)

        test_dic = self.dic_Cname2Hnames.copy()
        # Remove a random carbon
        del test_dic["C26"]
        assert not core.is_allHs_present(self.defop, self.dic_lipid, test_dic)

        test_dic2 = self.dic_Cname2Hnames.copy()
        # Change the number of hydrogens attached to a carbon
        test_dic2["CA2"] = ('HA21', 'HA22') # should be = ('HA21', 'HA22', 'HA23')
        assert not core.is_allHs_present(self.defop, self.dic_lipid, test_dic2)
