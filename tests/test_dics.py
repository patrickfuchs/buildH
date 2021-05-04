"""
Unit tests for buildH_calcOP.

Test functions from module init_dics
"""

import pathlib
import pytest

import MDAnalysis as mda

from buildh import lipids
from buildh import init_dics


DIR_DATA = "test_data"
PATH_ROOT_DATA = pathlib.Path(__file__).parent / DIR_DATA

# Ignore some MDAnalysis warnings for this test file
pytestmark = pytest.mark.filterwarnings('ignore::UserWarning')

# Called once before running the following tests
@pytest.fixture(scope='session')
def inputs():
    """Define some input parameters for the tests."""
    # Input parameters
    path_data = PATH_ROOT_DATA / "Berger_POPC"

    pdb = path_data / "10POPC.pdb"
    defop = path_data / "OP_def_BergerPOPC.def"
    lipids_tops = lipids.read_lipids_topH([lipids.PATH_JSON/"Berger_POPC.json"])
    dic_lipid = lipids_tops["Berger_POPC"]

    universe_woH = mda.Universe(str(pdb))
    dic_atname2genericname = init_dics.make_dic_atname2genericname(defop)
    return {'universe':universe_woH, 'defop':defop, 'dic_lipid':dic_lipid,
            'dic_atname2genericname': dic_atname2genericname}


def test_make_dic_atname2genericname(inputs):
    """Test for make_dic_atname2genericname().

    Parameters
    ----------
    inputs : function
        Callback which generates input data.
    """
    dic = init_dics.make_dic_atname2genericname(inputs['defop'])

    assert dic[('C1', 'H11')]   == 'gamma1_1'
    assert dic[('C32', 'H322')] == 'g1_2'
    assert dic[('C39', 'H392')] == 'palmitoyl_C5b'
    assert dic[('C17', 'H171')] == 'oleoyl_C2a'
    assert dic[('CA2', 'HA23')] == 'oleoyl_C18c'

    # Test non existent file
    with pytest.raises(ValueError) as err:
        init_dics.make_dic_atname2genericname("none")
    assert "Something went wrong when parsing the order parameter definition in file" in str(err.value)


def test_init_dic_OP(inputs):
    """Test for init_dic_OP().

    Parameters
    ----------
    inputs : function
        Callback which generates input data.
    """
    dic_OP, dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(inputs['universe'],
                                                                    inputs['dic_atname2genericname'],
                                                                    inputs['dic_lipid']['resname'])
    # Number of Order parmeters
    assert len(dic_OP) == 82
    for key in [('C1', 'H11'), ('C44', 'H442'), ('C17', 'H171'), ('CA1', 'HA11')]:
        assert key  in dic_OP.keys()
    # Number of lipid molecules
    assert len(dic_OP[('C37', 'H371')]) == 10

    # Number of lipid molecules
    assert len(dic_corresp_numres_index_dic_OP) == 10
    assert dic_corresp_numres_index_dic_OP[2] == 0
    assert dic_corresp_numres_index_dic_OP[11] == 9


def test_make_dic_Cname2Hnames(inputs):
    """Test for make_dic_Cname2Hnames().

    Parameters
    ----------
    inputs : function
        Callback which generates input data.
    """
    dic_OP, _ = init_dics.init_dic_OP(inputs['universe'], inputs['dic_atname2genericname'],
                                      inputs['dic_lipid']['resname'])

    dic = init_dics.make_dic_Cname2Hnames(dic_OP)

    # Number of Carbons with attached hydrogens
    assert len(dic) == 40
    assert dic['C1']  == ('H11', 'H12', 'H13')
    assert dic['C13'] == ('H131',)
    assert dic['C49'] == ('H491', 'H492')
    assert dic['CA2'] == ('HA21', 'HA22', 'HA23')
