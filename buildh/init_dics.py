"""Module to initialize dictionaries used in the program."""
import collections

# For debugging.
# TODO: Remove it after implement logging feature
DEBUG=False


def make_dic_atname2genericname(filename):
    """
    Make a dict of correspondance between generic H names and PDB names.

    This dict will look like the following: {('C1', 'H11'): 'gamma1_1', ...}.
    Useful for outputing OP with generic names (such as beta1, beta 2, etc.).
    Such files can be found on the NMRlipids MATCH repository:
    https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs.

    Parameters
    ----------
    filename : str
        Filename containing OP definition
        (e.g. `Berger_POPC.def`).

    Returns
    -------
    Ordered dictionary
        Keys are tuples of (C, H) name, values generic name (as described
        above in this docstring). The use of an ordered dictionary ensures
        we get always the same order in the output OP.
    """
    dic = collections.OrderedDict()
    try:
        with open(filename, "r") as f:
            for line in f:
                # TODO: This line might have to be changed if the file contains more than
                # 4 columns.
                name, _, C, H = line.split()
                dic[(C, H)] = name
    except:
        raise ValueError("Something went wrong when parsing the order parameter definition in "
                          "file {}".format(filename))
    return dic


def init_dic_OP(universe_woH, dic_atname2genericname, resname):
    """Initialize the dictionary of result (`dic_op`).

    Initialize also the dictionary of correspondance
    between residue number (resid) and its index in dic_OP (`dic_corresp_numres_index_dic_OP`).

    To calculate the error, we need to first average over the
    trajectory, then over residues.
    Thus in dic_OP, we want for each key a list of lists, for example:
    OrderedDict([
                 (('C1', 'H11'), [[], [], ..., [], []]),
                 (('C1', 'H12'), [[], ..., []]),
                 ...
                 ])
    Thus each sublist will contain OPs for one residue.
    e.g. ('C1', 'H11'), [[OP res 1 frame1, OP res1 frame2, ...],
                         [OP res 2 frame1, OP res2 frame2, ...], ...]

    Parameters
    ----------
    universe_woH : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    dic_atname2genericname: ordered dictionary
        dict of correspondance between generic H names and PDB names.
    resname: str
        The name of the lipid.

    Returns
    -------
    ordered dictionary
        Each key of this dict is a couple carbon/H, and at the beginning it
        contains an empty list.
    dictionary
        contains the correspondance between the residue number and
        its index in dic_op
    """
    dic_OP = collections.OrderedDict()


    # Get list of residue id from the lipid name
    all_resids = universe_woH.select_atoms( f"resname {resname}").residues.resids
    nb_residus = len(all_resids)

    # Each key contain a list which contains a number of list equals to
    # the number of residus
    for key in dic_atname2genericname:
        dic_OP[key] = [[] for _ in range(nb_residus)]

    # We also need the correspondance between residue number (resid) and
    # its index in dic_OP.
    # the index will always start at 0 and goes to the number of residus = range(nb_residus)
    dic_corresp_numres_index_dic_OP = dict(zip(all_resids, range(nb_residus)))

    if DEBUG:
        print("Initial dic_OP:", dic_OP)
        print("dic_corresp_numres_index_dic_OP:", dic_corresp_numres_index_dic_OP)

    return dic_OP, dic_corresp_numres_index_dic_OP


def make_dic_Cname2Hnames(dic_OP):
    """Initialize a dictionary of hydrogens bound to a carbon.

    Each key is a carbon and the value is a list of hydrogens bound to it.

    Parameters
    ----------
    dic_OP : ordered dictionary
        structure holding the calculated order parameter.

    Returns
    -------
    dictionary
        the constructed dictionary.
    """
    dic = {}
    for Cname, Hname in dic_OP.keys():
        if Cname not in dic:
            dic[Cname] = (Hname,)
        else:
            dic[Cname] += (Hname,)
    if DEBUG:
        print("dic_Cname2Hnames contains:", dic)
    return dic
