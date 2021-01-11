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
    "ForceField_lipidname" (e.g Berger_POPC, CHARMM_DMPC)

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
    UserWarning
        When a file doesn't have a correct format.
    """
    lipids_topH = {}
    for filename in filenames:
        filenam_path = pathlib.Path(filename)
        with open(filenam_path) as json_file:
            topol = json.load(json_file)
            # make sure at least 'resname' key exists
            if "resname" not in topol:
                raise UserWarning(f"{filenam_path} is in a bad format.")

            # Retrieve forcefield and lipid name from the filename
            ff, lipid_name = filenam_path.stem.split("_")

            # Generate keys by combining forcefield, the lipid name from the filename
            # and the possibles others lipids names from the json 'resname' attribut.
            # set() remove duplicates
            resnames = set([lipid_name] + topol["resname"])

            # Create distinct topology dictionaries for each resname supported in the json.
            for resname in resnames:
                data = topol.copy()
                data['resname'] = resname
                key = f"{ff}_{resname}"
                lipids_topH[key] = data

    return lipids_topH
