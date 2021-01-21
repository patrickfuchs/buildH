"""Provide functions to write the Order Parameters into files."""

import numpy as np

# For debugging.
# TODO: Remove it after implement logging feature
DEBUG=False

def pandasdf2pdb(df):
    """Return a string in PDB format from a pandas dataframe.

    Parameters
    ----------
    df : pandas dataframe with columns "atnum", "atname", "resname", "resnum",
         "x", "y", "z"

    Returns
    -------
    str
        A string representing the PDB.
    """
    s = ""
    chain = ""
    for _, row_atom in df.iterrows():
        atnum, atname, resname, resnum, x, y, z = row_atom
        atnum = int(atnum)
        resnum = int(resnum)
        # See for pdb format:
        # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html.
        # "alt" means alternate location indicator
        # "code" means code for insertions of residues
    	# "seg" means segment identifier
        # "elt" means element symbol
        if len(atname) == 4:
            s += ("{record_type:6s}{atnum:5d} {atname:<4s}{alt:1s}{resname:>4s}"
                  "{chain:1s}{resnum:>4d}{code:1s}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
                  "{occupancy:>6.2f}{temp_fact:>6.2f}          {seg:<2s}{elt:>2s}\n"
                  .format(record_type="ATOM", atnum=atnum, atname=atname, alt="",
                          resname=resname, chain=chain, resnum=resnum, code="",
                          x=x, y=y, z=z, occupancy=1.0, temp_fact=0.0, seg="",
                          elt=atname[0]))
        else:
            s += ("{record_type:6s}{atnum:5d}  {atname:<3s}{alt:1s}{resname:>4s}"
                  "{chain:1s}{resnum:>4d}{code:1s}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
                  "{occupancy:>6.2f}{temp_fact:>6.2f}          {seg:<2s}{elt:>2s}\n"
                  .format(record_type="ATOM", atnum=atnum, atname=atname, alt="",
                          resname=resname, chain=chain, resnum=resnum, code="",
                          x=x, y=y, z=z, occupancy=1.0, temp_fact=0.0, seg="",
                          elt=atname[0]))
    return s



def write_OP(fileout, dic_atname2genericname, dic_OP, resname):
    """Write the order parameters into a file.

    The output style comes from J. Melcr's script from NMRLipids project
    (https://github.com/NMRLipids/MATCH/blob/master/scripts/calcOrderParameters.py)


    Parameters
    ----------
    fileout: str
        name of the output file
    dic_atname2genericname: ordered dictionary
        dict of correspondance between generic H names and PDB names.
    dic_OP : ordered dictionary
        Each key of this dict is a couple carbon/H with the OP values as a list.
    resname : str
        lipid residue name taken from the json file.
    """
    with open(fileout, "w") as f:

        f.write("# {:18s} {:7s} {:5s} {:5s}  {:7s} {:7s} {:7s}\n"
                .format("OP_name", "resname", "atom1", "atom2", "OP_mean",
                        "OP_stddev", "OP_stem"))
        f.write("#-------------------------------"
                "-------------------------------------\n")
        # Loop over each pair (C, H).
        for Cname, Hname in dic_atname2genericname.keys():
            name = dic_atname2genericname[(Cname, Hname)]
            if DEBUG:
                print("Pair ({}, {}):".format(Cname, Hname))
            # Cast list of lists to a 2D-array. It should have dimensions
            # (nb_lipids, nb_frames).
            ### Thus each sublist will contain OPs for one residue.
            ### e.g. ('C1', 'H11'), [[OP res 1 frame1, OP res1 frame2, ...],
            ###                      [OP res 2 frame1, OP res2 frame2, ...],
            ####                     ...]
            a = np.array(dic_OP[(Cname, Hname)])
            if DEBUG:
                print("Final OP array has shape (nb_lipids, nb_frames):", a.shape)
                print()
            # General mean over lipids and over frames (for that (C, H) pair).
            mean = np.mean(a)
            # Average over frames for each (C, H) pair.  Because of how the
            # array is organized (see above), we need to average horizontally
            # (i.e. using axis=1).
            # means is a 1D-array with nb_lipids elements.
            means = np.mean(a, axis=1)
            # Calc standard deviation and STEM (std error of the mean).
            std_dev = np.std(means)
            stem = np.std(means) / np.sqrt(len(means))
            f.write("{:20s} {:7s} {:5s} {:5s} {: 2.5f} {: 2.5f} {: 2.5f}\n"
                    .format(name, resname, Cname, Hname, mean,
                            std_dev, stem))


def write_OP_alternate(fileout, universe_woH, dic_OP, resname):
    """Write the order parameters into a file with an alternate style.

    This style comes from A. Pineiro's script from NMRLipids project
    (https://github.com/NMRLipids/MATCH/blob/master/scratch/opAAUA_prod.py)


    Parameters
    ----------
    fileout: str
        name of the output file
    universe_woH : MDAnalysis universe instance
        This is the universe *without* hydrogen.
    dic_OP : ordered dictionary
        Each key of this dict is a couple carbon/H with the OP values as a list.
    resname : str
        lipid residue name taken from the json file.
    """
    with open(fileout, "w") as f:
        f.write("Atom_name  Hydrogen\tOP\t      STD\t   STDmean\n")
        list_unique_Cnames = []
        for Cname, Hname in dic_OP.keys():
            if Cname not in list_unique_Cnames:
                list_unique_Cnames.append(Cname)
        # Order of carbons is similar to that in the PDB.
        list_unique_Cnames_ordered = []
        selection = f"resname {resname}"
        for atom in universe_woH.select_atoms(selection).residues[0].atoms:
            if atom.name in list_unique_Cnames:
                list_unique_Cnames_ordered.append(atom.name)
        # Now write output.
        for Cname in list_unique_Cnames_ordered:
            cumulative_list_for_that_carbon = []
            for i, Hname in enumerate([H for C, H in dic_OP.keys() if C == Cname]):
                cumulative_list_for_that_carbon += dic_OP[Cname, Hname]
                a = np.array(dic_OP[Cname, Hname])
                mean = np.mean(a)
                means = np.mean(a, axis=1)
                std_dev = np.std(means)
                stem = np.std(means) / np.sqrt(len(means))
                if i == 0:
                    f.write("{:>7s}\t{:>8s}  {:10.5f}\t{:10.5f}\t{:10.5f}\n"
                            .format(Cname, "HR", mean, std_dev, stem))
                elif i == 1:
                    f.write("{:>7s}\t{:>8s}  {:10.5f}\t{:10.5f}\t{:10.5f}\n"
                            .format("", "HS", mean, std_dev, stem))
                elif i == 2:
                    f.write("{:>7s}\t{:>8s}  {:10.5f}\t{:10.5f}\t{:10.5f}\n"
                            .format("", "HT", mean, std_dev, stem))
            a = np.array(cumulative_list_for_that_carbon)
            mean = np.mean(a)
            means = np.mean(a, axis=1)
            std_dev = np.std(means)
            stem = np.std(means) / np.sqrt(len(means))
            f.write("{:>7s}\t{:>8s}  {:10.5f}\t{:10.5f}\t{:10.5f}\n\n"
                    .format("", "AVG", mean, std_dev, stem))
