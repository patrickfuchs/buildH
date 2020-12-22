#! /usr/bin/env python
# -*- coding: utf-8 -*-

import collections

# For debugging.
# TODO: Remove it after implement logging feature
DEBUG=False


def make_dic_atname2genericname(filename):
    """Make a dict of correspondance between generic H names and PDB names.

    This dict will look like the following: {('C1', 'H11'): 'gamma1_1', ...}.
    Useful for outputing OP with generic names (such as beta1, beta 2, etc.).
    Such files can be found on the NMRlipids MATCH repository:
    https://github.com/NMRLipids/MATCH/tree/master/scripts/orderParm_defs.

    Parameters
    ----------
    filename : str
        Filename containing OP definition
        (e.g. `order_parameter_definitions_MODEL_Berger_POPC.def`).

    Returns
    -------
    Ordered dictionnary
        Keys are tuples of (C, H) name, values generic name (as described
        above in this docstring). The use of an ordered dictionnary ensures
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
        raise UserWarning("Can't read order parameter definition in "
                          "file {}".format(filename))
    return dic

def init_dic_OP(universe_woH, dic_atname2genericname, dic_lipid):
    """TODO Complete docstring.
    """
    ### To calculate the error, we need to first average over the
    ### trajectory, then over residues.
    ### Thus in dic_OP, we want for each key a list of lists, for example:
    ### OrderedDict([
    ###              (('C1', 'H11'), [[], [], ..., [], []]),
    ###              (('C1', 'H12'), [[], ..., []]),
    ###              ...
    ###              ])
    ### Thus each sublist will contain OPs for one residue.
    ### e.g. ('C1', 'H11'), [[OP res 1 frame1, OP res1 frame2, ...],
    ###                      [OP res 2 frame1, OP res2 frame2, ...], ...]
    dic_OP = collections.OrderedDict()
    # We also need the correspondance between residue number (resnum) and
    # its index in dic_OP.
    dic_corresp_numres_index_dic_OP = {}
    # Create these sublists by looping over each lipid.
    for key in dic_atname2genericname:
        dic_OP[key] = []
        # Get lipid name.
        resname = dic_lipid["resname"]
        selection = "resname {}".format(resname)
        # Loop over each residue on which we want to calculate the OP on.
        for i, residue in enumerate(universe_woH.select_atoms(selection).residues):
            dic_OP[key].append([])
            dic_corresp_numres_index_dic_OP[residue.resid] = i
    if DEBUG:
        print("Initial dic_OP:", dic_OP)
        print("dic_corresp_numres_index_dic_OP:", dic_corresp_numres_index_dic_OP)
    return dic_OP, dic_corresp_numres_index_dic_OP


def make_dic_Cname2Hnames(dic_OP):
    """TODO Complete Docstring.
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
