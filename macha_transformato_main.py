#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 12:03:58 2022

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""

import sys
import os
import glob

from macha_transformato_functions import Preparation, CharmmManipulation

################################################################################
### VARIABLES/SETTINGS
################################################################################

parent_dir = "."
original_dir = "original"
ligands_dir = "ligands"
input_ext = "pdb"  # for testing - should be pdb
cgenff_path = "/site/raid2/johannes/programs/silcsbio/silcsbio.2022.1/cgenff/cgenff"

################################################################################
# MAIN (WORK)
################################################################################
# GET LIGAND ID(S) AND PUT THESE IN A LIST
# DOC: If called with an argument, a single ligand id will be processed.
# DOC: If no argument is given, the folder ligands_dir is searched for files
# DOC: with the extenstion input_ext and the names of these files will be
# DOC: used.

# Check for ligand as argument, and if found, run in single ligand mode
try:
    ligand_id = sys.argv[1]
    # Check for existence of ligands/"ligandid"
    if os.path.exists(f"{parent_dir}/{ligands_dir}/{ligand_id}.{input_ext}"):
        pass
    else:
        sys.exit(
            f"Input file: {parent_dir}/{ligands_dir}/{ligand_id}.{input_ext} not found!"
        )

# If no argument is given, run in multiple ligand mode
except IndexError:
    ligand_id = None  # None in particular

# Create the master list of the ligand/system names and fill it based on
# the findings above
ligand_ids = []
if ligand_id == None:
    for ifile in glob.glob(f"{parent_dir}/{original_dir}/*.{input_ext}"):
        ligand_ids.append(os.path.splitext(os.path.basename(ifile))[0])  #
elif ligand_id != "":
    ligand_ids.append(ligand_id)
else:
    sys.exit(f"Unknown ligand (file) name error with name: {ligand_id}")

###############################################################################

# ITERATE THROUGH LIGANDS
for ligand_id in ligand_ids:

    print(f"Processing ligand {ligand_id}")

    preparation = Preparation(
        parent_dir=parent_dir, ligand_id=ligand_id, original_dir=original_dir
    )

    # Make a Transformato style folder structure below a folder bearing
    # the name of the ligand
    preparation.makeTFFolderStructure()
    segids = preparation.createCRDfiles()

    # Get the toppar stream from a local CGenFF binary
    preparation.getTopparFromLocalCGenFF(cgenff_path=cgenff_path)
    preparation.copyREStocomplex()

    # COPY TEMPLATE FROM TEMPLATE FOLDER
    charmmManipulation = CharmmManipulation(
        parent_dir=parent_dir, ligand_id=ligand_id, original_dir=original_dir
    )

    charmmManipulation.manipulateToppar(preparation.resname)
    charmmManipulation.copyINPfiles()
    charmmManipulation.prepareStep1(segids, env = "waterbox")
    charmmManipulation.prepareStep1(segids, env = "complex")

    # MODIFY INPUT FILES FOR COMPLEX AND WATERBOX (FOR THIS LIGAND)

    # Suggested usage of inputFileInserter
    # - concatenate the blocks to be inserted using the respective functions
    # - then give the case ["stream toppar.str"]
    # - tell inputFileInserter to insert the block _after_ the case is found, so;
    # inputFileInserter(pathtoinputfile, ["stream toppar.str"], [concatenated_blocks], [False])
    # Note that each argument is a list with a single value - it needs to be this way since the function was made to make
    # several such insertions in the same file, at different places.
