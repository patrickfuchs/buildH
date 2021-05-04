"""
Unit tests for buildH.

Test functions from module UI
"""

import pathlib
import pytest


from buildh import UI, lipids, BuildHError
import sys

DIR_DATA = "test_data"
PATH_ROOT_DATA = pathlib.Path(__file__).parent / DIR_DATA

# Ignore some MDAnalysis warnings for this test file
pytestmark = pytest.mark.filterwarnings('ignore::UserWarning')

# path for the Berger POPC files
PATH_DATA = PATH_ROOT_DATA / "Berger_POPC"


class TestMain:
    """Test class for the main function of buildH."""

    # Arguments of the main function
    args = {
        "coord_file"        : str(PATH_DATA / "2POPC.pdb"),
        "traj_file"         : None,
        "def_file"          : str(PATH_DATA / "OP_def_BergerPOPC.def"),
        "out_file"          : "OP_buildH.out",
        "prefix_traj_ouput" : None,
        "begin"             : None,
        "end"               : None,
        "dic_lipid"         : None
    }

    # Default output used for assessement
    stdout_output = "Results written to OP_buildH.out"

    def setup_class(self):
        """Initialize attributes."""
        # Create correct lipid topology dict
        lipids_tops = lipids.read_lipids_topH([lipids.PATH_JSON/"Berger_POPC.json"])
        dic_lipid = lipids_tops["Berger_POPC"]
        self.args["dic_lipid"] = dic_lipid


    def test_main_minimal(self,  capsys):
        """Test main with minimal arguments."""
        UI.main(**self.args)
        captured = capsys.readouterr().out
        assert self.stdout_output in captured


    def test_main_output(self, capsys):
        """Test main with user defined output name."""
        args = self.args.copy()
        args["out_file"] = "text.txt"
        UI.main(**args)
        captured = capsys.readouterr().out
        assert "Results written to text.txt" in captured


    def test_main_subdef(self, capsys):
        """Test main with partial def file."""
        args = self.args.copy()
        args["def_file"] = str(PATH_DATA / "OP_def_HP_BergerPOPC.def")
        UI.main(**args)
        captured = capsys.readouterr().out
        assert self.stdout_output in captured


    def test_main_traj(self, capsys):
        """Test main with trajectory."""
        args = self.args.copy()
        args["traj_file"] = str(PATH_DATA / "2POPC.xtc")
        UI.main(**args)
        captured = capsys.readouterr().out
        assert self.stdout_output in captured
        assert "Dealing with frame 10 at 10000.0 ps." in captured


    def test_main_traj_slice(self, capsys):
        """Test main with sliced trajectory."""
        args = self.args.copy()
        args["traj_file"] = str(PATH_DATA / "2POPC.xtc")
        args["end"] = 3000
        UI.main(**args)
        captured = capsys.readouterr().out
        assert self.stdout_output in captured
        assert "Dealing with frame 3 at 3000.0 ps." in captured
        # Make sure we stop at frame 3
        assert "Dealing with frame 10 at 10000.0 ps." not in captured


    def test_main_traj_output(self, capsys):
        """Test main with trajectory and output trajectory."""
        args = self.args.copy()
        args["traj_file"] = str(PATH_DATA / "2POPC.xtc")
        args["prefix_traj_ouput"] = "test"
        UI.main(**args)
        captured = capsys.readouterr().out
        assert self.stdout_output in captured
        assert "Writing new pdb with hydrogens." in captured
        assert "Writing trajectory with hydrogens in xtc file." in captured


    def test_fail_main_coord_topology_mismatch(self):
        """ Test when coord file and topology doesn't match."""
        args = self.args.copy()

        lipids_tops = lipids.read_lipids_topH([lipids.PATH_JSON/"CHARMM_POPC.json"])
        dic_lipid = lipids_tops["CHARMM_POPC"]
        args["dic_lipid"] = dic_lipid

        with pytest.raises(BuildHError) as err:
            UI.main(**args)
        assert "The topology chosen does not match the structure provided" in str(err.value)


    def test_fail_main_coord_def_mismatch(self):
        """ Test when coord file and def file doesn't match."""
        args = self.args.copy()
        args["def_file"] = str(PATH_DATA / "op_wrong1.def")
        with pytest.raises(BuildHError) as err:
            UI.main(**args)
        assert f"Atoms defined in {args['def_file']} are missing in the structure" in str(err.value)


    def test_fail_main_subdef_traj(self,):
        """Test main with partial def file and a output trajectory. Must fail"""
        args = self.args.copy()
        args["def_file"] = str(PATH_DATA / "OP_def_HP_BergerPOPC.def")
        args["traj_file"] = str(PATH_DATA / "2POPC.xtc")
        args["prefix_traj_ouput"] = "test"
        with pytest.raises(BuildHError) as err:
            UI.main(**args)
        assert "Error on the number of H's to rebuild" in str(err.value)


class TestCLI:
    """Test class for the entry point of buildH.

    Mimic the command line arguments.
    """

    # Arguments of the CLI
    common_args = [
        "buildH",
        "-c", str(PATH_DATA / "2POPC.pdb"),
        "-d", str(PATH_DATA / "OP_def_BergerPOPC.def")
    ]

    def test_CLI_minimal(self, capsys):
        """Test working CLI with minimal arguments."""
        sys.argv = self.common_args + ["-l", "Berger_POPC"]
        UI.entry_point()
        captured = capsys.readouterr().out
        assert "Results written to OP_buildH.out" in captured


    def test_CLI_traj(self, capsys):
        """Test working CLI with all trajectory arguments."""
        sys.argv = (self.common_args + ["-t", str(PATH_DATA / "2POPC.xtc")]
                    + ["-l", "Berger_POPC"] + ["-o", "out.txt"]
                    + ["-opx", "base"] + ["-b", "1000", "-e", "10000"])
        UI.entry_point()
        captured = capsys.readouterr().out
        assert "Results written to out.txt" in captured
        assert "Dealing with frame 10 at 10000.0 ps." in captured
        assert "Writing new pdb with hydrogens." in captured
        assert "Writing trajectory with hydrogens in xtc file." in captured


    def test_CLI_user_json(self, capsys):
        """Test working CLI with JSON topology file."""
        sys.argv = (self.common_args + ["-l", "Berger_POPC"]
                    + ["-lt", str(PATH_ROOT_DATA / "Berger_POPC.json")])
        UI.entry_point()
        captured = capsys.readouterr().out
        assert "Results written to OP_buildH.out" in captured


    def test_fails_CLI_lipidtype(self, capsys):
        """Fail tests when passing a wrong lipid type."""
        sys.argv = self.common_args + ["-l", "PPHA"]
        with pytest.raises(SystemExit) as err:
            UI.entry_point()
        # Make sur the exception is thrown
        assert err.type == SystemExit
        assert "Lipid PPHA is not supported" in capsys.readouterr().err


    def test_fails_CLI_slice(self, capsys):
        """Fail tests when passing a slice option without a trajectory."""
        sys.argv =  sys.argv = self.common_args + ["-l", "Berger_POPC", "-e", "1000"]
        with pytest.raises(SystemExit) as err:
            UI.entry_point()
        assert err.type == SystemExit
        assert "Slicing is only possible with a trajectory file." in capsys.readouterr().err


class TestLaunch:
    """Test class for the launch function of buildH.

    This is the function called when using buildH as a module."""

    # Arguments of the main function
    args = {
        "coord_file"        : str(PATH_DATA / "2POPC.pdb"),
        "def_file"          : str(PATH_DATA / "OP_def_BergerPOPC.def"),
        "lipid_type"        : "Berger_POPC",
        "traj_file"         : None,
        "out_file"          : "OP_buildH.out",
        "prefix_traj_ouput" : None,
        "begin"             : None,
        "end"               : None,
        "lipid_jsons"       : None
    }

    # Default output used for assessement
    stdout_output = "Results written to OP_buildH.out"


    def test_launch_minimal(self, capsys):
        """Test launch with minimal arguments."""
        UI.launch(**self.args)
        captured = capsys.readouterr().out
        assert "Results written to OP_buildH.out" in captured

    def test_launch_traj(self, capsys):
        """Test launch with all trajectory arguments."""
        args = self.args.copy()
        args["traj_file"] = str(PATH_DATA / "2POPC.xtc")
        args["out_file"] = "out.txt"
        args["prefix_traj_ouput"] = "basename"
        args["begin"] = 0
        args["end"] = 10000
        UI.launch(**args)
        captured = capsys.readouterr().out
        assert "Results written to out.txt" in captured
        assert "Dealing with frame 10 at 10000.0 ps." in captured
        assert "Writing new pdb with hydrogens." in captured
        assert "Writing trajectory with hydrogens in xtc file." in captured


    def test_launch_user_json(self, capsys):
        """Test launch with JSON topology file."""
        args = self.args.copy()
        args["lipid_jsons"] = [str(PATH_ROOT_DATA / "Berger_POPC.json")]
        UI.launch(**args)
        captured = capsys.readouterr().out
        assert "Results written to OP_buildH.out" in captured


    def test_fail_launch_json(self):
        """Test fail launch with a wrong argument for JSON topology file."""
        args = self.args.copy()
        # Pass a string instead of a list
        args["lipid_jsons"] = str(PATH_ROOT_DATA / "Berger_POPC.json")
        with pytest.raises(TypeError) as err:
            UI.launch(**args)
        assert "a list is expected for argument 'lipid_jsons'" in str(err.value)


    def test_fail_launch_file(self):
        """Test fail launch with a non existant file."""
        args = self.args.copy()
        # Pass a string instead of a list
        args["traj_file"] = "nofile.xtc"
        with pytest.raises(FileNotFoundError) as err:
            UI.launch(**args)
        assert "nofile.xtc does not exist." in str(err.value)


    def test_fails_launch_lipidtype(self, capsys):
        """Fail tests when passing a wrong lipid type."""
        args = self.args.copy()
        args["lipid_type"] = "PPHA"
        with pytest.raises(BuildHError) as err:
            UI.launch(**args)
        assert err.type == BuildHError
        assert "Lipid PPHA is not supported" in str(err.value)
