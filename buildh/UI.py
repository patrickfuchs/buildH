"""User Interface module.

It contains functions used for the CLI mode
and the one for the 'module mode'.
"""

import argparse
import pickle
import pathlib
import sys

import numpy as np
import MDAnalysis as mda

from . import __version__
from . import lipids
from . import init_dics
from . import core
from . import writers
from . import utils


class BuildHError(Exception):
    """Raised when BuildH returns a error.

    This should be used when using BuildH as a module.
    """

    pass


def isfile(path):
    """Callback for checking file existence.

    This function checks if `path` is an existing file.
    If not, raise an error. Else, return the `path`.

    Parameters
    ----------
    path : str
        The path to be checked.

    Returns
    -------
    str
        The validated path.
    """
    source = pathlib.Path(path)
    if not pathlib.Path.is_file(source):
        if pathlib.Path.is_dir(source):
            msg = f"{source} is a directory."
        else:
            msg = f"{source} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return path


def parse_cli():
    """Handle the user parameters from the command line."""
    # Retrieve list of supported lipids
    lipids_files = [f for f in lipids.PATH_JSON.iterdir() if f.is_file()]
    lipids_tops = lipids.read_lipids_topH(lipids_files)
    lipids_supported_str = ", ".join(lipids_tops.keys())

    message = """This program builds hydrogens and calculates the order
    parameters (OP) from a united-atom trajectory of lipids. If -opx is
    requested, pdb and xtc output files with hydrogens are created but
    OP calculation will be slow. If no trajectory output is requested
    (no use of flag -opx), it uses a fast procedure to build hydrogens
    and calculate the OP.
    """
    epilog = f"The list of supported lipids (-l option) are: {lipids_supported_str}. "
    epilog += "More documentation can be found at https://buildh.readthedocs.io."
    parser = argparse.ArgumentParser(description=message,
                                     epilog=epilog)
    parser.add_argument('-v', '--version', action='version', version=__version__)
    # Avoid tpr for topology cause there's no .coord there!
    parser.add_argument("-c", "--coord", type=isfile, required=True,
                        help="Coordinate file (pdb or gro format).")
    parser.add_argument("-t", "--traj", type=isfile,
                        help="Input trajectory file. Could be in XTC, TRR or DCD format.")
    parser.add_argument("-l", "--lipid", type=str, required=True,
                        help="Combinaison of ForceField name and residue name "
                        "for the lipid to calculate the OP on (e.g. Berger_POPC)."
                        "It must match with the internal topology files or the one(s) supplied."
                        "A list of supported terms is printed when calling the help.")
    parser.add_argument("-lt", "--lipid_topology", type=isfile, nargs='+',
                        help="User topology lipid json file(s).")
    parser.add_argument("-d", "--defop", required=True, type=isfile,
                        help="Order parameter definition file. Can be found on "
                        "https://github.com/patrickfuchs/buildH/tree/master/def_files.")
    parser.add_argument("-opx", "--opdbxtc", help="Base name for trajectory "
                        "output with hydrogens. File extension will be "
                        "automatically added. For example -opx trajH will "
                        "generate trajH.pdb and trajH.xtc. "
                        "So far only xtc is supported.")
    parser.add_argument("-o", "--out", help="Output file name for storing "
                        "order parameters. Default name is OP_buildH.out.",
                        default="OP_buildH.out")
    parser.add_argument("-b", "--begin", type=int,
                        help="The first frame (ps) to read from the trajectory.")
    parser.add_argument("-e", "--end", type=int,
                        help="The last frame (ps) to read from the trajectory.")
    # parser.add_argument("-pi", "--pickle", type=str,
    #                     help="Output pickle filename. The structure pickled is a dictonnary "
    #                     "containing for each Order parameter, "
    #                     "the value of each lipid and each frame as a matrix")
    options = parser.parse_args()


    # Check coord file extension.
    if not options.coord.endswith("pdb") and not options.coord.endswith("gro"):
        parser.error("Coordinates file must be given in .pdb or .gro format.")

    # Use only the user lipid topologies
    if options.lipid_topology:
        lipids_files = [pathlib.Path(f) for f in options.lipid_topology]
        # Regenerate lipid topologies dictionary
        try:
            lipids_tops = lipids.read_lipids_topH(lipids_files)
        except ValueError as e:
            parser.error(e)
        # Regenerate str list of supported lipids.
        lipids_supported_str = ", ".join(lipids_tops.keys())

    # Check residue name validity.
    # Get the dictionary with helper info using residue name (options.lipid
    # argument).
    try:
        lipids_info = lipids_tops[options.lipid]
    except KeyError:
        parser.error(f"Lipid {options.lipid} is not supported. "
                     f"List of supported lipids are: {lipids_supported_str}")

    # Slicing only makes sense with a trajectory
    if not options.traj and (options.begin or options.end):
        parser.error("Slicing is only possible with a trajectory file.")

    return options, lipids_info


def entry_point():
    """Correspond to the entry point `buildH`."""
    # Parse arguments.
    args, dic_lipid = parse_cli()

    try:
        main(args.coord, args.traj, args.defop, args.out, args.opdbxtc, dic_lipid, args.begin, args.end)
    except BuildHError as e:
        sys.exit(e)


def launch(coord_file, def_file, lipid_type, traj_file=None, out_file="OP_buildH.out", prefix_traj_ouput=None, begin=None, end=None, lipid_jsons=None):
    """Launch BuildH calculations.

    This is the only function which can be called inside a Python script to use BuildH as a module.
    It checks the different arguments and call the main function.

    Parameters
    ----------
    coord_file : str
        Coordinate file. Only .pdb and .gro files are currently supported.
    def_file : str
        Order parameter definition file.
    lipid_type : str
        Combinaison of ForceField name and residue name for the lipid to calculate the OP on (e.g. Berger_POPC).
        It must match with the internal topology files or the one(s) supplied.
    traj_file : str, optional
        Trajectory file (could be in XTC, TRR or DCD format), by default None.
    out_file : str, optional
        Output file name for storing order parameters, by default "OP_buildH.out".
    prefix_traj_ouput : str, optional
        Base name for trajectory output with hydrogens.
        File extension will be automatically added.
        By default None.
    begin : int, optional
        The first frame (ps) to read from the trajectory, by default None.
    end : int, optional
        The last frame (ps) to read from the trajectory, by default None.
    lipid_jsons : list, optional
        User topology lipid json file(s), by default None.

    Raises
    ------
    FileNotFoundError
        When either coord_file, def_file or the traj_file is missing.
    TypeError
        When lipid_jsons is not a list.
    BuildHError
        When something went wront during calculation.
    """
    # Check files
    for file in [coord_file, traj_file, def_file]:
        if file is not None:
            source = pathlib.Path(file)
            if not pathlib.Path.is_file(source):
                raise FileNotFoundError(f"{source} does not exist.")

    # Construct available lipid topologies with only the user topology files.
    if lipid_jsons:
        if not isinstance(lipid_jsons, list):
            raise TypeError(f"a list is expected for argument 'lipid_jsons'")
        lipids_files = [pathlib.Path(f) for f in lipid_jsons]
        # Generate lipid topologies dictionary
        try:
            lipids_tops = lipids.read_lipids_topH(lipids_files)
        except ValueError as e:
            raise BuildHError(e)
    # Generate lipid topologies dictionary from internal topology files.
    else:
        lipids_files = [f for f in lipids.PATH_JSON.iterdir() if f.is_file()]
        lipids_tops = lipids.read_lipids_topH(lipids_files)

    try:
        dic_lipid = lipids_tops[lipid_type]
    except KeyError:
        raise BuildHError(f"Lipid {lipid_type} is not supported. "
                               f"List of supported lipids are: " + ", ".join(lipids_tops.keys()))


    try:
        main(coord_file, traj_file, def_file, out_file, prefix_traj_ouput, dic_lipid, begin, end)
    except BuildHError as e:
        raise e


def main(coord_file, traj_file, def_file, out_file, prefix_traj_ouput, dic_lipid, begin=None, end=None):
    """Main function of BuildH.

    It takes care of all the necessary steps to compute the Order Parameter :
      - create MDAnalysis Universe
      - create internal dictionaries
      - Check topology
      - Reconstruct the hydrogens
      - Compute OP and write the result in the output file.
      - If asked, write the trajectory with the hydrogens.

    Note
    ----
    This function shouldn't be called directly.
    Either it's called by the entry point or by the `launch` function.

    Parameters
    ----------
    coord_file : str
        Coordinate file. Only .pdb and .gro files are currently supported.
    traj_file : str
        Trajectory file (could be in XTC, TRR or DCD format). Can be None.
    def_file : str
        Order parameter definition file.
    out_file : str
        Output file name for storing order parameters.
    prefix_traj_ouput : str
        Base name for trajectory output with hydrogens.
        File extension will be automatically added.
        Can be None.
    dic_lipid : dict
        Lipid Topology for the reconstruction of the hydrogens.
    begin : int
        The first frame (ps) to read from the trajectory, by default None.
    end : int
        The last frame (ps) to read from the trajectory, by default None.

    Raises
    ------
    BuildHError
        When something went wront during calculation.
    """
    # Create universe without H.
    print("Constructing the system...")
    if traj_file:
        try:
            universe_woH = mda.Universe(coord_file, traj_file)
            begin, end = utils.check_slice_options(universe_woH, begin, end)
            traj_file = True
        except IndexError:
            raise BuildHError("Slicing options are not correct.")
        except:
            raise BuildHError(f"Can't create MDAnalysis universe with files {coord_file} and {traj_file}.")
    else:
        try:
            universe_woH = mda.Universe(coord_file)
            begin = 0
            end = 1
            traj_file = False
        except:
            raise BuildHError(f"Can't create MDAnalysis universe with file {coord_file}.")


    # Initialize dic for storing OP.
    # Init dic of correspondance : {('C1', 'H11'): 'gamma1_1',
    # {('C1', 'H11'): 'gamma1_1', ...}.
    try:
        dic_atname2genericname = init_dics.make_dic_atname2genericname(def_file)
    except ValueError as e:
        raise BuildHError(e)
    # Initialize dic_OP (see function init_dic_OP() for the format).
    dic_OP, dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(universe_woH,
                                                                    dic_atname2genericname,
                                                                    dic_lipid['resname'])
    # Initialize dic_Cname2Hnames.
    dic_Cname2Hnames = init_dics.make_dic_Cname2Hnames(dic_OP)


    # Now, a few checks have to be performed to ensure the lipid topology chosen (-l option),
    # the structure provided and the def file provider are coherent with each others.

    # Check if the lipid topology match the the structure.
    if not lipids.check_topology(universe_woH, dic_lipid):
        raise BuildHError(f"The topology chosen does not match the structure provided in {coord_file}")

    # Check if atoms name in the def file are present in the structure.
    atoms_name = [heavy_atom for (heavy_atom, _) in dic_atname2genericname.keys()]
    if not utils.check_def_file(universe_woH, dic_lipid['resname'], atoms_name):
        raise BuildHError(f"Atoms defined in {def_file} are missing in the structure {coord_file}.")

    # Check the def file and the topology are coherent.
    if not utils.check_def_topol_consistency(dic_Cname2Hnames, dic_lipid):
        raise BuildHError(f"Atoms defined in {def_file} are not consistent with the chosen topology.")


    print("System has {} atoms".format(len(universe_woH.coord)))

    # If traj output files are requested.
    # NOTE Here, we need to reconstruct all Hs. Thus the op definition file (passed
    #  with arg -d) needs to contain all possible C-H pairs !!!
    if prefix_traj_ouput:

        if utils.is_allHs_present(def_file, dic_lipid, dic_Cname2Hnames):
            core.gen_coordinates_calcOP(prefix_traj_ouput, universe_woH, dic_OP, dic_lipid,
                                        dic_Cname2Hnames, dic_corresp_numres_index_dic_OP,
                                        begin, end, traj_file)
        else:
            raise BuildHError(f"Error on the number of H's to rebuild. An output trajectory has been "
                     f"requestest but {def_file} doesn't contain all hydrogens to rebuild.")

    # If no traj output file requested, use fast indexing to speed up OP
    # calculation. The function fast_build_all_Hs() returns nothing, dic_OP
    # is modified in place.
    if not prefix_traj_ouput:
        core.fast_build_all_Hs_calc_OP(universe_woH, begin, end, dic_OP,
                                       dic_lipid, dic_Cname2Hnames)


    # Output to a file.
    writers.write_OP(out_file, dic_atname2genericname,
                            dic_OP, dic_lipid['resname'])
    print(f"Results written to {out_file}")

    # # Pickle results
    # if args.pickle:
    #     with open(args.pickle, "wb") as f:
    #         # Pickle the dic using the highest protocol available.
    #         pickle.dump(dic_OP, f, pickle.HIGHEST_PROTOCOL)
    #         print("dictionary pickled and written to {}".format(args.pickle))
    #     #  To unpickle
    #     #with open("OP.pickle", "rb") as f:
    #     #    dic_OP = pickle.load(f)
