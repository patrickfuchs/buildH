"""
Module for the lipid topology json files.

This module contains functions for parsing the json files.
"""

import pathlib
import json

# Directory name of the json files
JSON_DIR = "lipids"
# Absolute path of the json files
PATH_JSON = pathlib.Path(__file__).parent / JSON_DIR

def read_lipids_topH(filenames):
    """Generate a list of lipid hydrogen topologies.

    This function read a list of json files containing
    the topology of a united-atom lipid for reconstructing
    the missing hydrogens.

    The list topologies is stored as dictionary were the key is
    "ForceField_lipidname" (e.g Berger_POPC, CHARMM_DMPC).
    If there is multiple lipid residue name in the key 'resname' of the json,
    it will create one dictionary for one element of this key.

    Parameters
    ----------
    filenames : list of str
        List of json files containing topologies.

    Returns
    -------
    dictionary of dictionary
        the lipid topologies.

    Raises
    ------
    ValueError
        When a file doesn't have a correct format.
    """
    lipids_tops = {}
    for filename in filenames:
        filenam_path = pathlib.Path(filename)
        with open(filenam_path) as json_file:
            try:
                topol = json.load(json_file)
            except Exception as e:
                raise ValueError(f"{filenam_path} is in a bad format: " + str(e)) from e
            # make sure at least 'resname' key exists
            if "resname" not in topol:
                raise ValueError(f"{filenam_path} is in a bad format: keyword 'resname' is missing.")

            # Retrieve forcefield and lipid name from the filename
            try:
                ff, lipid_name = filenam_path.stem.split("_")
            except ValueError as e:
                raise ValueError(f"{filenam_path} has an incorrect name. "
                                 "It should be Forcefield_LipidName.json") from e

            # Generate keys by combining forcefield, the lipid name from the filename
            # and the possibles others lipids names from the json 'resname' attribut.
            # set() remove duplicates
            resnames = set([lipid_name] + topol["resname"])

            # Create distinct topology dictionaries for each resname supported in the json.
            for resname in resnames:
                data = topol.copy()
                data['resname'] = resname
                key = f"{ff}_{resname}"
                lipids_tops[key] = data

    return lipids_tops


def check_topology(universe, lipid_top):
    """Check if the topology `lipid_top` is coherent with the structure in `universe`.

    This function check first the lipid residue name then the atoms in the topology.
    This function return false if there is one missing in the structure.
    Print also an error message.

    Parameters
    ----------
    universe : MDAnalysis universe instance
    lipid_top : dict
        lipid topology for hydrogen.

    Returns
    -------
    Boolean
        whether the structure and the topology match
    """
    # Check first the lipid residue name
    resname = lipid_top['resname']
    lipid_atoms = universe.select_atoms( f"resname {resname}")
    if len(lipid_atoms) == 0:
        print(f"No lipid '{resname}' found in the topology.")
        return False

    # Remove first key 'resname'
    carbon_atoms = list(lipid_top.keys())[1::]
    #retrieve all atom names in the system
    all_names = set(universe.select_atoms(f"resname {resname}").names)
    if not set(carbon_atoms).issubset(all_names):
        miss_atoms = ",".join(set(carbon_atoms) - all_names)
        print(f"Some atoms ({miss_atoms}) from topology are not found in your system.")
        return False

    return True
