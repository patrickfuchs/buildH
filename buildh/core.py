"""Module holding the core functions."""

import pandas as pd
import MDAnalysis as mda
import MDAnalysis.coordinates.XTC as XTC

from . import hydrogens
from . import geometry as geo
from . import writers


# For debugging.
# TODO: Remove it after logging feature implementation.
DEBUG=False


def buildHs_on_1C(atom, H_type, helper1, helper2, helper3=None):
    """Build 1, 2 or 3 H on a given carbon.

    This function is a wrapper which takes the coordinates of the helpers
    and call the function that builds 1, 2 or 3 H.

    Parameters
    ----------
    atom : numpy 1D-array
        Central atom on which we want to reconstruct the hydrogen.
    H_type: str
        The type of H to build. It could be 'CH2', 'CH', 'CHdoublebond' or 'CH3'
        see dic_lipds.py
    helper1 : numpy 1D-array
        First neighbor of central atom.
    helper2 : numpy 1D-array
        Second neighbor of central atom.
    helper3 : numpy 1D-array
        Third neighbor of central atom.

    Returns
    -------
    tuple of numpy 1D-arrays
        Each element of the tuple is a numpy 1D-array containing 1, 2 or 3
        reconstructed hydrogen(s).
        !!! IMPORTANT !!! This function *should* return a tuple even if
        there's only one H that has been rebuilt.
    """
    if H_type == "CH2":
        H1_coor, H2_coor = hydrogens.get_CH2(atom, helper1, helper2)
        return (H1_coor, H2_coor)
    elif H_type == "CH":
        # If we reconstruct a single H, we have a 3rd helper.
        #helper3_coor = sel("name {0}".format(helper3_name))[0].position
        H1_coor = hydrogens.get_CH(atom, helper1, helper2,
                                   helper3)
        return (H1_coor,)
    elif H_type == "CHdoublebond":
        H1_coor = hydrogens.get_CH_double_bond(atom, helper1,
                                               helper2)
        return (H1_coor,)
    elif H_type == "CH3":
        H1_coor, H2_coor, H3_coor = hydrogens.get_CH3(atom,
                                                      helper1, helper2)
        return (H1_coor, H2_coor, H3_coor)
    else:
        raise UserWarning("Wrong code for typeofH2build, expected 'CH2', 'CH'"
                          ", 'CHdoublebond' or 'CH3', got {}."
                          .format(H_type))


def build_system_hydrogens(universe_woH, dic_lipid, dic_Cname2Hnames, dic_lipid_indexes):
    """Build a new system *with* hydrogens.

    The function takes the MDAnalysis universe *without* hydrogens, reconstructs all hydrogens
    and returns a pandas dataframe. This latter will be used later to build a new
    MDAnalysis universe with H.

    Notes
    -----
    There is no simple way to create a new MDAnalysis universe directly.

    Parameters
    ----------
    universe_woH : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    dic_lipid : dictionary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.
    dic_Cname2Hnames : dictionary
        This dict gives the correspondance Cname -> Hname. It is a dict of
        tuples. If there is more than 1 H for a given C, they need to be
        *ordered* like in the PDB. e.g. for CHARMM POPC :
        {'C13': ('H13A', 'H13B', 'H13C'), ..., 'C33': ('H3X', 'H3Y'),
          ..., 'C216': ('H16R', 'H16S'), ...}
    dic_lipids_with_indexes : dictionary
        The dictionary made in function make_dic_lipids_with_indexes().

    Returns
    -------
    pandas dataframe
        contains the system *with* hydrogens.
    """
    # The list newrows will be used to store the new molecule *with* H.
    newrows = []
    # Counter for numbering the new mlcs with H.
    new_atom_num = 1
    # Loop over all atoms in the universe without H.
    for atom in universe_woH.atoms:
        # Make resnum resets to 0 when it exceeds 9999.
        resnum = atom.resnum-((atom.resnum//10000)*10000)
        resname = atom.resname
        name = atom.name
        # Append atom to the new list.
        # 0      1       2        3       4  5  6
        # atnum, atname, resname, resnum, x, y, z
        newrows.append([new_atom_num, name, resname, resnum]
                        + list(atom.position))
        new_atom_num += 1
        # Build new H(s)?
        if (atom.name in dic_lipid and atom.residue.resname == dic_lipid["resname"]):
            # Retrieve helpers coordinates
            # helperX_ix is the index of the helper inside one residue.
            if len(dic_lipid_indexes[atom.name]) == 6:
                typeofH2build, _, _, _, helper1_ix, helper2_ix = dic_lipid_indexes[atom.name]
                helper3_coor = None
            else:
                typeofH2build, _, _, _, _, helper1_ix, helper2_ix, helper3_ix = dic_lipid_indexes[atom.name]
                helper3_coor = atom.residue.atoms[helper3_ix].position

            helper1_coor = atom.residue.atoms[helper1_ix].position
            helper2_coor = atom.residue.atoms[helper2_ix].position

            # Build Hs and store them in a list of numpy 1D-arrays Hs_coor.
            # The "s" in Hs_coor means there can be more than 1 H:
            # For CH2, Hs_coor will contain: [H1_coor, H2_coor].
            # For CH3, Hs_coor will contain: [H1_coor, H2_coor, H3_coor].
            # For CH, Hs_coor will contain: [H1_coor].
            # For CHdoublebond, Hs_coor will contain: [H1_coor].
            Hs_coor = buildHs_on_1C(atom.position, typeofH2build,
                                    helper1_coor, helper2_coor, helper3_coor)

            # Loop over Hs_coor (H_coor is a 1D-array with the 3 coors of 1 H).
            for i, H_coor in enumerate(Hs_coor):
                # Retrieve name of newly built H.
                Hname = dic_Cname2Hnames[atom.name][i]
                # Add them to newrows.
                newrows.append([new_atom_num, Hname, resname, resnum]
                                + list(H_coor))
                new_atom_num += 1

    # Create a dataframe to store the mlc with added hydrogens.
    new_df_atoms = pd.DataFrame(newrows, columns=["atnum", "atname",
                                                  "resname", "resnum",
                                                  "x", "y", "z"])
    return new_df_atoms


###
### The next function build_all_Hs_calc_OP())
### build new H, calculate the order parameter and write the new traj with Hs
### to an output file (e.g. .xtc, etc).
### Note: it is slow, it shouldn't be used if the user doesn't want to
###       write the trajectory. Instead, fast_build_all_Hs() should be used.
###
def build_all_Hs_calc_OP(universe_woH, ts, dic_lipid, dic_Cname2Hnames, universe_wH, dic_OP,
                         dic_corresp_numres_index_dic_OP, dic_lipid_indexes):
    """Build all hydrogens and calculates order parameters for one frame.

    This function loop overs *all* atoms of the universe_woH in order to update
    the atom coordinates and the new H built into the universe_wH.

    The function also calculates the order parameter.
    The coordinates of the universe *with* H are updated in place.
    The order parameter is also added in place (within dic_OP dictionary).

    Notes
    -----
    This function is slow, thus it shall be used when one wants
    to create a trajectory with H (such as .xtc or whatever format).

    This function assumes all possible C-H pairs are present in the .def
    file (with -d option). They are needed since we want to build an xtc with
    the whole system. If one is interested in calculating only a subset of OPs,
    please use the function fast_build_all_Hs_calc_OP() instead.

    Parameters
    ----------
    universe_woH : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    ts : Timestep instance
        the current timestep with the coordinates
    dic_lipid : dictionary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.
    dic_Cname2Hnames : dictionary
        This dict gives the correspondance Cname -> Hname. It is a dict of
        tuples. If there is more than 1 H for a given C, they need to be
        *ordered* like in the PDB. e.g. for CHARMM POPC :
        {'C13': ('H13A', 'H13B', 'H13C'), ..., 'C33': ('H3X', 'H3Y'),
          ..., 'C216': ('H16R', 'H16S'), ...}
    universe_wH : MDAnalysis universe instance (optional)
        This is the universe *with* hydrogens.
    dic_OP : ordered dictionary
        Each key of this dict is a couple carbon/H, and at the beginning it
        contains an empty list, e.g.
        OrderedDict([ ('C1', 'H11): [], ('C1', 'H12'): [], ... ])
        See function init_dic_OP() below to see how it is organized.
    dic_corresp_numres_index_dic_OP : dictionary
        This dict should contain the correspondance between the numres and
        the corresponding index in dic_OP. For example {..., 15: 14, ...} means
        the residue numbered 15 in the PDB has an index of 14 in dic_OP.
    """
    # We will need the index in the numpy array for updating coordinates
    # in the universe with H.
    row_index_coor_array = 0
    resid = -9999

    # Loop over all atoms in the universe without H.
    for atom in universe_woH.atoms:
        # Update the position of the current atom in the universe with H.
        universe_wH.coord.positions[row_index_coor_array, :] = atom.position
        row_index_coor_array += 1

        # Build new H(s)?
        if (atom.name in dic_lipid and atom.residue.resname == dic_lipid["resname"]):

            # Retrieve the index of the first atom in the current residue
            # Test to avoid refreshing it at every step of the loop
            if resid != atom.residue.resid:
                resid = atom.residue.resid
                ix_first_atom_res = atom.residue.atoms[0].ix

            # Retrieve helpers coordinates
            if len(dic_lipid_indexes[atom.name]) == 6:
                typeofH2build, _, _, _, helper1_ix, helper2_ix = dic_lipid_indexes[atom.name]
                helper3_coor = None
            else:
                typeofH2build, _, _, _, _, helper1_ix, helper2_ix, helper3_ix = dic_lipid_indexes[atom.name]
                helper3_coor = ts[helper3_ix + ix_first_atom_res]

            # Faster to retrieve the coordinates from ts than from universe_woH.atoms.positions
            helper1_coor = ts[helper1_ix + ix_first_atom_res]
            helper2_coor = ts[helper2_ix + ix_first_atom_res]

            # Build Hs and store them in a list of numpy 1D-arrays Hs_coor.
            # The "s" in Hs_coor means there can be more than 1 H:
            # For CH2, Hs_coor will contain: [H1_coor, H2_coor].
            # For CH3, Hs_coor will contain: [H1_coor, H2_coor, H3_coor].
            # For CH, Hs_coor will contain: [H1_coor].
            # For CHdoublebond, Hs_coor will contain: [H1_coor].
            Hs_coor = buildHs_on_1C(atom.position, typeofH2build,
                                    helper1_coor, helper2_coor, helper3_coor)

            # Loop over Hs_coor (H_coor is a 1D-array with the 3 coors of 1 H).
            for i, H_coor in enumerate(Hs_coor):
                # Retrieve name of newly built H.
                Hname = dic_Cname2Hnames[atom.name][i]
                ####
                #### We calculate here the order param on the fly :-D !
                ####
                if (atom.name, Hname) in dic_OP:
                    op = geo.calc_OP(atom.position, H_coor)
                    # We should get here the index of the residue in dic_OP.
                    # For that we can use dic_corresp_numres_index_dic_OP
                    # (key: resnum in pdb, value: index residue in dic_OP).
                    lipid_ix = dic_corresp_numres_index_dic_OP[atom.resid]
                    # OLD way: dic_OP[(atom.name, Hname)].append(op)
                    if (atom.name, Hname) in dic_OP:
                        dic_OP[(atom.name, Hname)][lipid_ix].append(op)
                    if DEBUG:
                        print(atom.name, H_coor, "OP:", op)

                # Update the position of the current H in the universe with H.
                universe_wH.coord.positions[row_index_coor_array, :] = H_coor
                row_index_coor_array += 1

            if dic_OP and DEBUG:
                print()
                print()
    if dic_OP and DEBUG:
        print("Final dic_OP:", dic_OP)
        print()


###
### The next 3 functions (get_indexes(), make_dic_lipids_with_indexes()
### and fast_build_all_Hs_calc_OP()) should be used when the
### user doesn't want an output trajectory.
### By using fast indexing to individual Catoms and helpers, they
### are much faster.
###
def get_indexes(atom, dic_lipid):
    """Return the index of helpers for a given carbon.

    Parameters
    ----------
    atom : MDAnalysis Atom instance
        This is an Atom instance of a carbon on which we want to build Hs.
    dic_lipid : dictionary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.

    Returns
    -------
    tuple of 2 or 3 int
        The tuple contains the index of the 2 (or 3) helpers for the atom that
        was passed as argument. (e.g. for atom C37 with index 99, the function
        returns a tuple containing 98 (index of C36 = helper 1) and 100 (index
        of C38=helper2).
    """
    # Get nb of H to build and helper names (we can have 2 or 3 helpers).
    if len(dic_lipid[atom.name]) == 3:
        typeofH2build, helper1_name, helper2_name = dic_lipid[atom.name]
    else:
        typeofH2build, helper1_name, helper2_name, helper3_name = dic_lipid[atom.name]
    # Get helper coordinates using atom, which an instance from Atom class.
    # atom.residue.atoms is a list of atoms we can select with
    # method .select_atoms().
    # To avoid too long line, we shorten its name to `sel`.
    sel = atom.residue.atoms.select_atoms
    helper1_ix = sel("name {}".format(helper1_name))[0].ix
    helper2_ix = sel("name {}".format(helper2_name))[0].ix
    if typeofH2build == "CH":
        # If we reconstruct a single H, we have a 3rd helper.
        helper3_ix = sel("name {0}".format(helper3_name))[0].ix
        return (helper1_ix, helper2_ix, helper3_ix)
    else:
        return (helper1_ix, helper2_ix)


def make_dic_lipids_with_indexes(universe_woH, dic_lipid, dic_OP):
    """Expand dic_lipid and adds the index of each atom and helper.

    IMPORTANT: the index of each atom/helper is given with respect to the
               first atom in that residue.
    For example, if we have a POPC where C1 is the first atom, and C50 the
    last one, we want in the end:
    {'C1': ('CH3', 'N4', 'C5', 0, 3, 4), ...,
     'C50': ('CH3', 'C49', 'C48', 49, 48, 47)}
    Where the 3 last int are the index (ix) of the atom, helper1, helper2
    (possibly helper3) with respect to the first atom.
    Thus for C1 : 0 is index of C1, N4 is 3 atoms away from C1 and C5 is 4
    atoms away from C1.
    For C50: C50 is 49 atoms away from C1, C49 is 48 atoms away from C1,
    C48 is 47 atoms away from C1.

    Parameters
    ----------
    universe_woH : MDAnalysis Universe instance
        The universe without hydrogens.
    dic_lipid : dictionary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.
    dic_OP : ordered dictionary
        Each key of this dict is a couple carbon/H, and at the beginning it
        contains an empty list, e.g.
        OrderedDict([ ('C1', 'H11): [], ('C1', 'H12'): [], ... ])
        See function init_dic_OP() below to see how it is organized.

    Returns
    -------
    dictionary
        The returned dictionary as described above in this docstring.
    """
    # Get lipid name.
    resname = dic_lipid["resname"]
    # Get resnum of the 1st lipid encountered in the system whose name
    # is `resname`.
    selection = "resname {}".format(resname)
    first_lipid_residue = universe_woH.select_atoms(selection).residues[0]
    resnum_1st_lipid = first_lipid_residue.resnum
    # Get name of 1st atom of that lipid.
    first_atom_name = first_lipid_residue.atoms[0].name
    # Get index of this atom.
    first_atom_ix = first_lipid_residue.atoms[0].ix
    if DEBUG:
        print("resname: {}, first encountered residue: {},\n"
              "resnum_1st_lipid: {}, first_atom_name: {}, first_atom_ix: {}"
              .format(resname, first_lipid_residue, resnum_1st_lipid,
                      first_atom_name, first_atom_ix))
        print()
    # Keep only carbons on which we want to build Hs.
    carbons2keep = []
    for Cname, _ in dic_OP:
        if Cname not in carbons2keep:
            carbons2keep.append(Cname)
    dic_lipids_with_indexes = {}
    for Cname in dic_lipid.keys():
        if Cname in carbons2keep:
            dic_lipids_with_indexes[Cname] = dic_lipid[Cname].copy()
    # Now add the helper indexes.
    # The reasonning is over one residue (e.g. POPC). We want to add (to the
    # dict) the index (ix) of each helper of a given carbon with respect to
    # the index of the first atom in that lipid residue.
    # Loop over each carbon on which we want to reconstruct Hs.
    for Cname in dic_lipids_with_indexes:
        # Loop over residues for a given Cname atom.
        selection = "resid {} and name {}".format(resnum_1st_lipid, Cname)
        for Catom in universe_woH.select_atoms(selection):
            # Get the (absolute) index of helpers.
            if dic_lipid[Cname][0] == "CH":
                helper1_ix, helper2_ix, helper3_ix = get_indexes(Catom, dic_lipid)
            else:
                helper1_ix, helper2_ix = get_indexes(Catom, dic_lipid)
            # If the first lipid doesn't start at residue 1 we must
            # substract the index of the first atom of that lipid.
            Catom_ix_inres = Catom.ix - first_atom_ix
            helper1_ix_inres = helper1_ix - first_atom_ix
            helper2_ix_inres = helper2_ix - first_atom_ix
            # Then add these indexes to dic_lipids_with_indexes.
            if dic_lipid[Cname][0] == "CH":
                helper3_ix_inres = helper3_ix - first_atom_ix
                tmp_tuple = (Catom_ix_inres, helper1_ix_inres,
                             helper2_ix_inres, helper3_ix_inres)
                dic_lipids_with_indexes[Cname] += tmp_tuple
            else:
                tmp_tuple = (Catom_ix_inres, helper1_ix_inres,
                             helper2_ix_inres)
                dic_lipids_with_indexes[Cname] += tmp_tuple
    if DEBUG:
        print("Everything is based on the following dic_lipids_with_indexes\n{}"
              .format(dic_lipids_with_indexes))
        print()
    return dic_lipids_with_indexes


def fast_build_all_Hs_calc_OP(universe_woH, begin, end,
                              dic_OP, dic_lipid, dic_Cname2Hnames):
    """Build Hs and calc OP using fast indexing.

    This function uses fast indexing to carbon atoms and helper atoms. It
    should be used when the user doesn't want any output traj with hydrogens.

    Parameters
    ----------
    universe_woH : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    begin: int
        index of the first frame of trajectory
    end: int
        index of the last frame of trajectory
    dic_OP : ordered dictionary
        Each key of this dict is a couple carbon/H, and at the beginning it
        contains an empty list, e.g.
        OrderedDict([ ('C1', 'H11): [], ('C1', 'H12'): [], ... ])
        See function init_dic_OP() below to see how it is organized.
    dic_lipid : dictionary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.
    dic_Cname2Hnames : dictionary
        This dict gives the correspondance Cname -> Hname. It is a dict of
        tuples. If there is more than 1 H for a given C, they need to be
        *ordered* like in the PDB. e.g. for CHARMM POPC :
        {'C13': ('H13A', 'H13B', 'H13C'), ..., 'C33': ('H3X', 'H3Y'),
          ..., 'C216': ('H16R', 'H16S'), ...}

    Returns
    -------
    None
        This function returns nothing, dic_OP is changed *in place*.
    """
    ###
    ### 1) Expand dic_lipids and store there helpers' index.
    ###
    ### We want {'C1': ('CH3', 'N4', 'C5', 0, 3, 4), ...,
    ###          'C50': ('CH3', 'C49', 'C48', 49, 48, 47)}
    ### Where the 3 last int are the index (ix) of the atom, helper1, helper2
    ### (possibly helper3) with respect to the first atom
    ### (e.g. 0 is index of C1, N4 is 3 atoms away from C1, etc).
    ###
    dic_lipids_with_indexes = make_dic_lipids_with_indexes(universe_woH,
                                                           dic_lipid, dic_OP)
    # Get lipid name.
    resname = dic_lipid["resname"]
    # Select first residue of that lipid.
    selection = "resname {}".format(resname)
    first_lipid_residue = universe_woH.select_atoms(selection).residues[0]
    # Get name of 1st atom of that lipid.
    first_atom_name = first_lipid_residue.atoms[0].name
    ###
    ### 2) Now loop over the traj, residues and Catoms.
    ### At each iteration build Hs and calc OP.
    ###
    # Loop over frames (ts is a Timestep instance).
    for ts in universe_woH.trajectory[begin:end]:
        if len(universe_woH.trajectory) == 1:
            print(f"Dealing with current frame.")
        else:
            print(f"Dealing with frame {ts.frame} at "
                  f"{universe_woH.trajectory.time} ps.")
        if DEBUG:
            print("Looping now over residues...")
            print()
        # Loop over the 1st atom of each lipid, which is equiv to loop *over
        # residues* (first_lipid_atom is an Atom instance, lipid_ix is an int
        # that will be used for storing OPs in dic_OP).
        selection = "resname {} and name {}".format(resname, first_atom_name)
        for lipid_ix, first_lipid_atom in enumerate(universe_woH.select_atoms(selection)):
            if DEBUG:
                print("Dealing with Cname", first_lipid_atom)
                print("which is part of residue", first_lipid_atom.residue)
                print("Now looping over atoms of this residue")
                print()
            # Get the index of this first atom.
            ix_first_atom_res = first_lipid_atom.ix
            # Now loop over each carbon on which we want to build Hs
            # (Cname is a string).
            for Cname in dic_lipids_with_indexes:
                # Get Cname and helpers coords.
                if len(dic_lipids_with_indexes[Cname]) == 6:
                    typeofH2build, _, _, Cname_ix, helper1_ix, helper2_ix = dic_lipids_with_indexes[Cname]
                    helper3_coor = None
                else:
                    typeofH2build, _, _, _, Cname_ix, helper1_ix, helper2_ix, helper3_ix = dic_lipids_with_indexes[Cname]
                    helper3_coor = ts[helper3_ix+ix_first_atom_res]
                helper1_coor = ts[helper1_ix+ix_first_atom_res]
                helper2_coor = ts[helper2_ix+ix_first_atom_res]
                Cname_position = ts[Cname_ix+ix_first_atom_res]
                if DEBUG:
                    print("Dealing with Cname", Cname)
                    sel = first_lipid_atom.residue.atoms.select_atoms
                    Cname_atom = sel("name {}".format(Cname))[0]
                    print(Cname_atom, Cname_atom.position)
                    if len(dic_lipid[Cname]) == 3:
                        _, helper1_name, helper2_name = dic_lipid[Cname]
                    else:
                        _, helper1_name, helper2_name, helper3_name = dic_lipid[Cname]
                    helper1_atom = sel("name {}".format(helper1_name))[0]
                    print("helper1", helper1_atom, helper1_atom.position)
                    helper2_atom = sel("name {}".format(helper2_name))[0]
                    print("helper2", helper2_atom, helper2_atom.position)
                    if len(dic_lipid[Cname]) == 4:
                        helper3_atom = sel("name {}".format(helper3_name))[0]
                        print("helper3", helper3_atom, helper3_atom.position)
                # Get newly built H(s) on that atom.
                Hs_coor = buildHs_on_1C(Cname_position, typeofH2build,
                                        helper1_coor, helper2_coor, helper3_coor)
                if DEBUG:
                    print("Cname_position with fast indexing:", Cname_position)
                    print("helper1_position with fast indexing:",
                          ts[helper1_ix+ix_first_atom_res])
                    print("helper2_position with fast indexing:",
                          ts[helper2_ix+ix_first_atom_res])
                    if len(dic_lipid[Cname]) == 4:
                        print("helper3_position with fast indexing:",
                              ts[helper3_ix+ix_first_atom_res])
                # To retrieve Hname, we need a counter.
                counter4Hname = 0
                # Loop over all Hs.
                for H_coor in Hs_coor:
                    # Retrieve name of newly built H.
                    Hname = dic_Cname2Hnames[Cname][counter4Hname]
                    # Calc and store OP for that couple C-H.
                    Cname_position = ts[Cname_ix+ix_first_atom_res]
                    op = geo.calc_OP(Cname_position, H_coor)
                    # Old way: dic_OP[(Cname, Hname)].append(op)
                    if (Cname, Hname) in dic_OP:
                        dic_OP[(Cname, Hname)][lipid_ix].append(op)
                    if DEBUG:
                        print(Hname, H_coor, "OP:", op)
                    # Increment counter4Hname for retrieving next H.
                    counter4Hname += 1
                if DEBUG:
                    print()
                    print()
    if DEBUG:
        print("Final dic_OP:", dic_OP)
        print()


def gen_coordinates_calcOP(basename, universe_woH, dic_OP, dic_lipid,
                           dic_Cname2Hnames, dic_corresp_numres_index_dic_OP,
                           begin, end, traj_file):
    """Generate coordinates files (pdb and/or xtc) with computed hydrogens
    and compute the order parameter.

    If `traj_file` is set to False, only a pdb file will be written.
    This depends whether or not the user supplied a trajectory file
    in the first place.

    Parameters
    ----------
    basename : str
        basename for the output coordinate file(s).
    universe_woH : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    dic_OP : ordered dictionary
        Each key of this dict is a couple carbon/H, and at the beginning it
        contains an empty list, e.g.
        OrderedDict([ ('C1', 'H11): [], ('C1', 'H12'): [], ... ])
        See function init_dic_OP() below to see how it is organized.
    dic_lipid : dictionary
        Comes from dic_lipids.py. Contains carbon names and helper names needed
        for reconstructing hydrogens.
    dic_Cname2Hnames : dictionary
        This dict gives the correspondance Cname -> Hname. It is a dict of
        tuples. If there is more than 1 H for a given C, they need to be
        *ordered* like in the PDB. e.g. for CHARMM POPC :
        {'C13': ('H13A', 'H13B', 'H13C'), ..., 'C33': ('H3X', 'H3Y'),
          ..., 'C216': ('H16R', 'H16S'), ...}
    dic_corresp_numres_index_dic_OP : dictionary
        This dict should contain the correspondance between the numres and
        the corresponding index in dic_OP.
    begin: int
        index of the first frame of trajectory
    end: int
        index of the last frame of trajectory
    traj_file : bool
        a trajectory output file has to be generated?
    """
    dic_lipids_with_indexes = make_dic_lipids_with_indexes(universe_woH, dic_lipid,
                                                           dic_OP)

    # Create filenames.
    pdbout_filename = basename + ".pdb"
    # Build a new universe with H.
    # Build a pandas df with H.
    new_df_atoms = build_system_hydrogens(universe_woH, dic_lipid, dic_Cname2Hnames,
                                          dic_lipids_with_indexes)
    # Create a new universe with H using that df.
    print("Writing new pdb with hydrogens.")
    # Write pdb with H to disk.
    with open(pdbout_filename, "w") as f:
        f.write(writers.pandasdf2pdb(new_df_atoms))
    # Then create the universe with H from that pdb.
    if len(universe_woH.trajectory) > 1:
        universe_wH = mda.Universe(pdbout_filename, dt=universe_woH.trajectory.dt)
    else:
        universe_wH = mda.Universe(pdbout_filename)

    #Do we need to generate a trajectory file ?
    if traj_file:
        xtcout_filename = basename + ".xtc"
        # Create an xtc writer.
        print("Writing trajectory with hydrogens in xtc file.")
        newxtc = XTC.XTCWriter(xtcout_filename, len(universe_wH.atoms))

        # 4) Loop over all frames of the traj *without* H, build Hs and
        # calc OP (ts is a Timestep instance).
        for ts in universe_woH.trajectory[begin:end]:
            if len(universe_woH.trajectory) == 1:
                print(f"Dealing with current frame.")
            else:
                print(f"Dealing with frame {ts.frame} at "
                      f"{universe_woH.trajectory.time} ps.")
            # Build H and update their positions in the universe *with* H (in place).
            # Calculate OPs on the fly while building Hs  (dic_OP changed in place).
            build_all_Hs_calc_OP(universe_woH, ts, dic_lipid, dic_Cname2Hnames,
                                universe_wH, dic_OP, dic_corresp_numres_index_dic_OP,
                                dic_lipids_with_indexes)
            # Update the box values into the new Universe
            universe_wH.trajectory[0].dimensions = ts.dimensions
            # Write new frame to xtc.
            newxtc.write(universe_wH)
        # Close xtc.
        newxtc.close()
    # if not, just compute OP in the fast way.
    else:
        fast_build_all_Hs_calc_OP(universe_woH, begin, end, dic_OP, dic_lipid, dic_Cname2Hnames)
