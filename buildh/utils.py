"""Module for various control functions."""

import numpy as np

def check_slice_options(system, first_frame=None, last_frame=None):
    """Verify the slicing options and return a range of frames in MDAnalysis.

    This function check whether the first frame and the last frame are consistent
    within themselves (`first_frame` cant be superior to `last_frame`) and
    with the trajectory supplied (if the trajectory starts at 1000ps, `first_frame`
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
    IndexError
        When the slicing options are out of range.
    """
    # From the trajectory, get the time of the first and last frame
    traj_first_frame = int(system.trajectory.time)
    traj_last_frame = int(system.trajectory.time +
                          system.trajectory.dt * (system.trajectory.n_frames - 1))

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


def check_atom(universe, res_name, atom_name):
    """Check if `atom_name` from residue `res_name` is present in `universe`.

    Parameters
    ----------
    universe : MDAnalysis universe instance
    res_name : str
        residue name
    atom_name : str
        atom name

    Returns
    -------
    Bool
        True if the atom is present. False otherwise.
    """
    if len(universe.select_atoms(f"resname {res_name} and name {atom_name}")) == 0:
        return False
    return True


def check_def_file(universe, res_name, atoms_name):
    """Check if atoms from the definition file are present in the structure in `universe`.

    This function return false if there is one missing in the structure.
    Print also an error message.

    Parameters
    ----------
    universe : MDAnalysis universe instance
    res_name : str
        lipid residue name
    atoms_name : list of str
        list of atom names

    Returns
    -------
    Bool
        True if all atoms are found in the structure. False otherwise.
    """
    # get all atom names of the res_name in the system
    all_names = set(universe.select_atoms(f"resname {res_name}").names)
    if not set(atoms_name).issubset(all_names):
        miss_atoms = ",".join(set(atoms_name) - all_names)
        print(f"Some atoms ({miss_atoms}) of residue {res_name} from definition "
                "file are not found in your system.")
        return False

    return True


def check_def_topol_consistency(dic_Cname2Hnames, lipid_top):
    """Check the consistency between the lipid topology and the def file.

    Ensure that the carbons in the def file are present in the topology.
    Ensure that all hydrogens of a given carbon as described in the
    topology are present in the def file.

    Parameters
    ----------
    dic_Cname2Hnames : dict
        This dict gives the correspondance Cname -> Hname. It is a dict of
        tuples extracted from the def file.
    lipid_top : dict
        lipid topology for hydrogen.

    Returns
    -------
    Bool
        True is it's coherent. False otherwise.
    """
    # Check if carbons in dic_Cname2Hnames keys are all present in the lipid
    # topology.
    if not set(dic_Cname2Hnames.keys()).issubset(lipid_top.keys()):
        miss_atoms = ",".join(set(dic_Cname2Hnames.keys()) - set(lipid_top.keys()))
        print(f"Some carbons ({miss_atoms}) from the definition file are not"
              "present in the json file.")
        return False

    # For each carbon in topology, make sure all hydrogens attached
    # are in the def file
    nb_Hs_expected = {'CH3': 3, 'CH2': 2, 'CH': 1, 'CHdoublebond': 1}
    for carbon, values in lipid_top.items():
        if carbon != "resname":
            H_type = values[0]
            nb_Hs_topol = nb_Hs_expected[H_type]
            if carbon in dic_Cname2Hnames: # Handle partial def file
                nb_Hs_def = len(dic_Cname2Hnames[carbon])
                if nb_Hs_def != nb_Hs_topol:
                    print(f"Carbon {carbon} from the definition file should contains "
                        f"{nb_Hs_topol} hydrogen(s), found {nb_Hs_def}.")
                    return False
    return True


def is_allHs_present(def_file, lipids_name, dic_ref_CHnames):
    """Check if all H's to be rebuild are present in the def file.

    Parameters
    ----------
    def_file : str
        Filename containing OP definition
        (e.g. `Berger_POPC.def`).
    lipids_name : dictionary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.
    dic_ref_CHnames: dictionary
        Contains all CH molecules

    Returns
    -------
    boolean
        whether or not all Hs are present.
    """
    # check that dic_ref_CHnames contains all possible C-H pairs.
    # NOTE The user has to take care that .def file has the right atom names !!!
    for atname in lipids_name.keys():
        if atname != "resname":
            # Check if carbon is present in the definition file.
            if atname not in dic_ref_CHnames:
                print("Error: When -opx option is used, the order param "
                      "definition file (passed with -d arg) must contain "
                      "all possible carbons on which we want to rebuild "
                      "hydrogens.")
                print("Found:", list(dic_ref_CHnames.keys()))
                print("Needs:", list(lipids_name.keys()))
                return False
            # Check that the 3 Hs are in there for that C.
            nbHs_in_def_file = len(dic_ref_CHnames[atname])
            tmp_dic = {"CH": 1, "CHdoublebond": 1, "CH2": 2, "CH3": 3}
            correct_nb_of_Hs = tmp_dic[lipids_name[atname][0]]
            if  correct_nb_of_Hs != nbHs_in_def_file:
                print("Error: When -opx option is used, the order param "
                      "definition file (passed with -d arg) must contain "
                      "all possible C-H pairs to rebuild.")
                print("Expected {} hydrogen(s) to rebuild for carbon {}, "
                      "got {} in definition file {}."
                      .format(correct_nb_of_Hs, atname,
                              dic_ref_CHnames[atname], def_file))
                return False

    return True
