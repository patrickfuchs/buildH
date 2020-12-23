"""Command line entry point for buildH."""

import argparse
import pickle
import pathlib

import numpy as np
import MDAnalysis as mda

from . import dic_lipids
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


def parse_cli():
    """
    Handle the user parameters from the command line.
    """

    message = """This program builds hydrogens and calculate the order
    parameters
    (OP) from a united-atom trajectory. If -opx is requested, pdb and xtc
    output files with hydrogens are created but OP calculation will be slow.
    If no trajectory output is requested (no use of flag -opx), it uses a
    fast procedure to build hydrogens and calculate the OP.
    """
    parser = argparse.ArgumentParser(description=message)
    # Avoid tpr for topology cause there's no .coord there!
    parser.add_argument("topfile", type=isfile,
                        help="Topology file (pdb or gro).")
    parser.add_argument("-x", "--xtc", type=isfile,
                        help="Input trajectory file in xtc format.")
    parser.add_argument("-l", "--lipid", type=str, required=True,
                        help="Residue name of lipid to calculate the OP on (e.g. POPC).")
    parser.add_argument("-d", "--defop", required=True, type=isfile,
                        help="Order parameter definition file. Can be found on "
                        "NMRlipids MATCH repository:"
                        "https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs")
    parser.add_argument("-opx", "--opdbxtc", help="Base name for trajectory "
                        "output with hydrogens. File extension will be "
                        "automatically added. For example -opx trajH will "
                        "generate trajH.pdb and trajH.xtc. "
                        "So far only xtc is supported.")
    parser.add_argument("-o", "--out", help="Output base name for storing "
                        "order parameters. Extention \".out\" will be "
                        "automatically added. Default name is OP_buildH.out.",
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
    # Check residue name validity.
    # Get the dictionnary with helper info using residue name (options.lipid
    # argument). Beware, this dict is then called `dic_lipid` *without s*,
    # while `dic_lipids.py` (with an s) is a module with many different dicts
    # (of different lipids) the user can choose.
    try:
        lipids_info = getattr(dic_lipids, options.lipid)
    except AttributeError:
        parser.error("Lipid dictionnary {} doesn't exist in dic_lipids.py".format(options.lipid))

    # Slicing only makes sense with a trajectory
    if not options.xtc and (options.begin or options.end):
        parser.error("Slicing is only possible with a trajectory file.")

    return options, lipids_info


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
        except:
            raise UserWarning("Can't create MDAnalysis universe with file {}"
                              .format(args.topfile))
    print("System has {} atoms".format(len(universe_woH.coord)))

    # 2) Initialize dic for storing OP.
    # Init dic of correspondance : {('C1', 'H11'): 'gamma1_1',
    # {('C1', 'H11'): 'gamma1_1', ...}.
    dic_atname2genericname = init_dics.make_dic_atname2genericname(args.defop)
    # Initialize dic_OP (see function init_dic_OP() for the format).
    dic_OP, dic_corresp_numres_index_dic_OP = init_dics.init_dic_OP(universe_woH,
                                                                    dic_atname2genericname,
                                                                    dic_lipid)
    # Initialize dic_Cname2Hnames.
    dic_Cname2Hnames = init_dics.make_dic_Cname2Hnames(dic_OP)

    # If traj output files are requested.
    # NOTE Here, we need to reconstruct all Hs. Thus the op definition file (passed
    #  with arg -d) needs to contain all possible C-H pairs !!!
    if args.opdbxtc:

        if core.is_allHs_present(args.defop, dic_lipid, dic_Cname2Hnames):
            core.gen_XTC_calcOP(args.opdbxtc, universe_woH, dic_OP, dic_lipid,
                                dic_Cname2Hnames, dic_corresp_numres_index_dic_OP,
                                begin, end)
        else:
            raise UserWarning("Error on the number of H's to rebuild.")

    # 6) If no traj output file requested, use fast indexing to speed up OP
    # calculation. The function fast_build_all_Hs() returns nothing, dic_OP
    # is modified in place.
    if not args.opdbxtc:
        core.fast_build_all_Hs_calc_OP(universe_woH, begin, end, dic_OP, dic_lipid, dic_Cname2Hnames)


    # Output to a file.
    writers.write_OP_jmelcr("{}.jmelcr_style.out".format(args.out), dic_atname2genericname,
                            dic_OP, dic_lipid)
    writers.write_OP_apineiro("{}.apineiro_style.out".format(args.out), universe_woH,
                              dic_OP, dic_lipid)
    print("Results written to {}".format(args.out))

    # Pickle results
    if args.pickle:
        with open(args.pickle, "wb") as f:
            # Pickle the dic using the highest protocol available.
            pickle.dump(dic_OP, f, pickle.HIGHEST_PROTOCOL)
            print("Dictionnary pickled and written to {}".format(args.pickle))
        #  To unpickle
        #with open("OP.pickle", "rb") as f:
        #    dic_OP = pickle.load(f)
