#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 12:03:58 2022

@author: Johannes Karwanoupoulos and Åsmund Kaupang
"""

import sys
import os
import glob

from macha_transformato_functions import *

################################################################################
### VARIABLES/SETTINGS
################################################################################

original_dir = "original"
ligands_dir = "ligands"
input_ext = "pdb"


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
    if os.path.exists(f"{ligands_dir}/{ligand_id}.{input_ext}"):
        pass
    else:
        sys.exit(f"Input file: {ligands_dir}/{ligand_id}.{input_ext} not found!")
        
# If no argument is given, run in multiple ligand mode
except IndexError:
    ligand_id = None #None in particular

# Create the master list of the ligand/system names and fill it based on
# the findings above
ligand_ids = []
if ligand_id == None:
    for ifile in glob.glob(f"{ligands_dir}/*.{input_ext}"):
        ligand_ids.append(#
                          os.path.splitext(os.path.basename(ifile))[0]
                          )
elif ligand_id != '':
    ligand_ids.append(ligand_id)
else:
    sys.exit(f"Unknown ligand (file) name error with name: {ligand_id}")
    
###############################################################################

    
# based on ligands make transformato directories
# run cgenff on ligands
# copy from template dir to each
# modify input step1

