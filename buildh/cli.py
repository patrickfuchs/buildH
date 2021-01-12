"""Command line entry point for buildH."""

import argparse
import pickle
import pathlib
import sys

import numpy as np
import MDAnalysis as mda

from . import lipids
from . import init_dics
from . import core
from . import writers

# For debugging.
DEBUG = False
# For pickling results (useful for future analyses, e.g. drawing distributions).
PICKLE = False

def isfile(path):
    """Callback for checking file existence.

    This function checks if path is an existing file.
    If not, raise an error. Else, return the path.

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


def check_slice_options(system, first_frame=None, last_frame=None):
    """Verify the slicing options given by the user and translate
    to it to a range of frame in MDAnalysis.

    This function check whether the first frame and the last frame are consistent
    within themselves (``first_frame`` cant be superior to ``last_frame``) and
    with the trajectory supplied (if the trajectory starts at 1000ps, ``first_frame``
    cant be equal to 0 for example).
    Then, the function translate the range from picosecond-time to the number of frame
    in MDanalysis.


    Parameters
    ----------
    system : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    first_frame : int
        the first frame to read (in ps)
    last_frame : int
        the last frame to read (in ps)

    Return
    ------
    tuple of int
        The number of first and last frame

    Raises
    ------
    """
    # From the trajectory, get the time of the first and last frame
    traj_first_frame = int(system.trajectory.time)
    traj_last_frame = int(system.trajectory.time + system.trajectory.dt * (system.trajectory.n_frames - 1))

    # If no bound is given, take the full trajectory
    if not first_frame and not last_frame:
        return (0, system.trajectory.n_frames)

    # If only one bound is given
    if not first_frame:
        first_frame = traj_first_frame
    if not last_frame:
        last_frame = traj_last_frame


    # Check abnormal range
    if first_frame < 0 or last_frame < 0:
        raise IndexError("Incorrect slice options.")
    if first_frame > last_frame:
        raise  IndexError("Incorrect slice options")

    # Check if the range fits into the range of the trajectory
    if first_frame < traj_first_frame or last_frame < traj_first_frame:
        raise  IndexError("Incorrect slice options")
    if first_frame > traj_last_frame or last_frame > traj_last_frame:
        raise  IndexError("Incorrect slice options")

    # Translate the time range into a number range.
    # Find the index of element in the list of frames (in ps) which has the minimum distance
    # from the first or last frame (in ps) given.
    frames = np.arange(traj_first_frame, traj_last_frame + 1, int(system.trajectory.dt))
    number_first_frame = (np.abs(frames - first_frame)).argmin()
    number_last_frame  = (np.abs(frames -  last_frame)).argmin()
    # Include last frame into account for slicing by adding 1
    number_last_frame = number_last_frame + 1

    return (number_first_frame, number_last_frame)

def check_lipid_present(universe_woH, lipid_name):
    """Check if the `lipid_name` residue name is present in `universe_woH`.

    Parameters
    ----------
    universe_woH : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    lipid_name : str
        lipid residue name.

    Returns
    -------
    Boolean
        whether or not it's present.
    """
    lipid_atoms = universe_woH.select_atoms( f"resname {lipid_name}")

    if len(lipid_atoms) == 0:
        return False
    return True

def parse_cli():
    """
    Handle the user parameters from the command line.
    """

    # Retrieve list of supported lipids
    lipids_files = [f for f in lipids.PATH_JSON.iterdir() if f.is_file()]
    lipids_tops = lipids.read_lipids_topH(lipids_files)
    lipids_supported_str = ", ".join(lipids_tops.keys())

    message = """This program builds hydrogens and calculate the order
    parameters
    (OP) from a united-atom trajectory. If -opx is requested, pdb and xtc
    output files with hydrogens are created but OP calculation will be slow.
    If no trajectory output is requested (no use of flag -opx), it uses a
    fast procedure to build hydrogens and calculate the OP.
    """
    epilog = f"The list of supported lipids (-l option) are: {lipids_supported_str}."
    parser = argparse.ArgumentParser(description=message,
                                     epilog=epilog)
    # Avoid tpr for topology cause there's no .coord there!
    parser.add_argument("topfile", type=isfile,
                        help="Topology file (pdb or gro).")
    parser.add_argument("-x", "--xtc", type=isfile,
                        help="Input trajectory file in xtc format.")
    parser.add_argument("-l", "--lipid", type=str, required=True,
                        help="Residue name of lipid to calculate the OP on (e.g. POPC).")
    parser.add_argument("-tl", "--lipid_topology", type=isfile, nargs='+',
                        help="User topology lipid json file(s). Mandatory to build hydrogens.")
    parser.add_argument("-d", "--defop", required=True, type=isfile,
                        help="Order parameter definition file. Can be found on "
                        "NMRlipids MATCH repository:"
                        "https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs")
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
    parser.add_argument("-pi", "--pickle", type=str,
                        help="Output pickle filename. The structure pickled is a dictonnary "
                        "containing for each Order parameter, the value of each lipid and each frame as a matric")
    options = parser.parse_args()

    # Top file is "options.topfile", xtc file is "options.xtc", pdb output file is
    # "options.pdbout", xtc output file is "options.xtcout".
    # Check topology file extension.
    if not options.topfile.endswith("pdb") and not options.topfile.endswith("gro"):
        parser.error("Topology must be given in pdb or gro format")

    # Append user lipid topologies to the default ones
    if (options.lipid_topology):
        lipids_files += [pathlib.Path(f) for f in options.lipid_topology]
        # Regenerate lipid topologies dictionary
        lipids_tops = lipids.read_lipids_topH(lipids_files)
        # Regenerate str list of supported lipids.
        lipids_supported_str = ", ".join(lipids_tops.keys())

    # Check residue name validity.
    # Get the dictionnary with helper info using residue name (options.lipid
    # argument).
    try:
        lipids_info = lipids_tops[options.lipid]
    except KeyError:
        parser.error(f"Lipid {options.lipid} is not supported. List of supported lipids are: {lipids_supported_str}")

    # Slicing only makes sense with a trajectory
    if not options.xtc and (options.begin or options.end):
        parser.error("Slicing is only possible with a trajectory file.")

    return options, lipids_info


def check_atoms_def(universe_woH, res_name, atoms_name):
    print(atoms_name)
    for atom_name in atoms_name:
        if not check_atom_def(universe_woH, res_name, atom_name):
            print(f"Atom {atom_name} from definition file is not found in your system.")
            return False

    return True


def check_atom_def(universe_woH, res_name, atom_name):

    if len(universe_woH.select_atoms(f"resname {res_name} and name {atom_name}")) == 0:
        return False
    return True


def main():
    """buildH main function."""

    # 1) Parse arguments.
    args, dic_lipid = parse_cli()

    # 2) Create universe without H.
    print("Constructing the system...")
    if args.xtc:
        try:
            universe_woH = mda.Universe(args.topfile, args.xtc)
            begin, end = check_slice_options(universe_woH, args.begin, args.end)
            traj_file = True
        except IndexError:
            raise UserWarning("Slicing options are not correct.") from None
        except:
            raise UserWarning("Can't create MDAnalysis universe with files {} "
                              "and {}".format(args.topfile, args.xtc)) from None
    else:
        try:
            universe_woH = mda.Universe(args.topfile)
            begin = 0
            end = 1
            traj_file = False
        except:
            raise UserWarning("Can't create MDAnalysis universe with file {}"
                              .format(args.topfile))


    # 2) Initialize dic for storing OP.
    # Init dic of correspondance : {('C1', 'H11'): 'gamma1_1',
    # {('C1', 'H11'): 'gamma1_1', ...}.
    dic_atname2genericname = init_dics.make_dic_atname2genericname(args.defop)
    # Initialize dic_OP (see function init_dic_OP() for the format).
    dic_OP, dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(universe_woH,
                                                                    dic_atname2genericname,
                                                                    dic_lipid['resname'])
    # Initialize dic_Cname2Hnames.
    dic_Cname2Hnames = init_dics.make_dic_Cname2Hnames(dic_OP)


    # Now, a few checks have to be performed to ensure the lipid topology chosen (-l option),
    # the structure provided and the def file provider are coherent with each others.

    # Check if the lipid name from topology (`lipids_info`) is present in the structure.
    if not check_lipid_present(universe_woH, dic_lipid['resname']):
        sys.exit(f"No lipid '{dic_lipid['resname']}' found in {args.topfile}.")

    # Check if the topology chosen is coherent with the structure.

    # Check if atoms names in the def file are present in the structure.
    atoms_def = [heavy_atom for (heavy_atom, _) in dic_atname2genericname.keys()]
    if not check_atoms_def(universe_woH, dic_lipid['resname'],atoms_def):
        sys.exit(f"Atoms defined in {args.defop} are missing in the structure {args.topfile}.")



    print("System has {} atoms".format(len(universe_woH.coord)))

    # If traj output files are requested.
    # NOTE Here, we need to reconstruct all Hs. Thus the op definition file (passed
    #  with arg -d) needs to contain all possible C-H pairs !!!
    if args.opdbxtc:

        if core.is_allHs_present(args.defop, dic_lipid, dic_Cname2Hnames):
            core.gen_coordinates_calcOP(args.opdbxtc, universe_woH, dic_OP, dic_lipid,
                                        dic_Cname2Hnames, dic_corresp_numres_index_dic_OP,
                                        begin, end, traj_file)
        else:
            raise UserWarning("Error on the number of H's to rebuild.")

    # 6) If no traj output file requested, use fast indexing to speed up OP
    # calculation. The function fast_build_all_Hs() returns nothing, dic_OP
    # is modified in place.
    if not args.opdbxtc:
        core.fast_build_all_Hs_calc_OP(universe_woH, begin, end, dic_OP, dic_lipid, dic_Cname2Hnames)


    # Output to a file.
    writers.write_OP(args.out, dic_atname2genericname,
                            dic_OP, dic_lipid['resname'])
    print(f"Results written to {args.out}")

    # Pickle results
    if args.pickle:
        with open(args.pickle, "wb") as f:
            # Pickle the dic using the highest protocol available.
            pickle.dump(dic_OP, f, pickle.HIGHEST_PROTOCOL)
            print("Dictionnary pickled and written to {}".format(args.pickle))
        #  To unpickle
        #with open("OP.pickle", "rb") as f:
        #    dic_OP = pickle.load(f)
