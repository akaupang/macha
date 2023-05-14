#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated in May 2023

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""
import sys
import argparse
from macha.functions import checkInput, Preparation
from macha.charmm_factory import CharmmManipulation
import logging

################################################################################
### VARIABLES/SETTINGS
################################################################################

parent_dir = "data"
original_dir = "original"
input_ext = "pdb"         # exclusive support for PDB files for now
protein_name = "protein"  # -> protein.pdb
cgenff_path = ""          # MUST BE SET BY USER

################################################################################
### LOGGING
################################################################################
# set up logging to file
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-15s %(levelname)-8s %(message)s',
                    datefmt='%d.%m.%Y %H:%M',
                    filename=f"{parent_dir}/macha_run.log",
                    filemode='a') # change to "w" if you want the log file to be 
                                  # overwritten instead of appended to

# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)

# set a format which is simpler for console use
#formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
#formatter = logging.Formatter('%(levelname)-1s %(message)s')
formatter = logging.Formatter('%(message)s')

# tell the handler to use this format
console.setFormatter(formatter)

# add the handler to the root logger
logging.getLogger().addHandler(console)

# # But we prefer module-specific logging
# preplog = logging.getLogger('Preparation')
# cmanlog = logging.getLogger('charmmFactory')

################################################################################
# MAIN (WORK)
################################################################################
# Introductions are in order
logging.info(f"{str(' '):>90}")
logging.info(f"======================================================================   m   a   c   h   a")
logging.info(f"{str(' '):>90}")
       
# ARGUMENT HANDLING
# Initialize
parser = argparse.ArgumentParser(description='Options')
# Add expected/possible arguments
parser.add_argument('-nc', '--nocomplex', action='store_true', help="Toggle production of water boxes only (no complexes)")
parser.add_argument('-nw', '--nowaterbox', action='store_true', help="Toggle production of complexes only (no water boxes)")
parser.add_argument('-rna', '--rna', action='store_true', help="Toggle handling of RNA")
# Parse arguments
args = parser.parse_args()

# Determine which environments/run types will be used
if args.nocomplex == True:
    #no_complex = True
    envs = ["waterbox"]
    logging.info("Attention: Only water boxes will be produced!")
    logging.info("")
elif args.nowaterbox == True:
    #no_waterbox = True
    envs = ["complex"]
    logging.info("Attention: Only complexes will be produced!")
    logging.info("")
else:
    #no_complex = False
    envs = ["waterbox", "complex"]

   
rna = False 
if args.rna == True:
    rna = True
    logging.info("Attention: Input type set to RNA. Other input types will produce errors!")
    logging.info("")
else:
    rna = False
    
# Check for CGenFF (except for RNA?)
if ((cgenff_path == "") and (rna == False)):
    logging.critical(   f"\tThis package requires cgenff for parameterization.\n"\
                        f"\tPlease install it in the active environment or point the routine\n"\
                        f"\tto the right path using the key cgenff_path='/path/to/cgenff' .\n"\
                        f"\t"
                    )
    sys.exit(1)

# Check for protein.pdb and if none is found, assume complexes are provided.
# If protein.pdb is found, assume ligands alone.
protein_id, ligand_ids = checkInput(
    parent_dir = parent_dir, 
    original_dir = original_dir, 
    protein_name = protein_name
    )

# ITERATE THROUGH LIGANDS
for ligand_id in ligand_ids:
    # Pretty printing
    logging.info(f"__________________________________________________________________________________________")
    logging.info(f"{str(' '):>90}")
    logging.info(f"Processing {ligand_id}")

    for env in envs:
        # Announce the current environment
        logging.info(f"_________________________________________________________________________________ {env.upper():>8}")
        logging.info(f"{str(' '):>90}")


        ########################################################################
        # Preparation of ligands
        # Instantiate class
        preparation = Preparation(
            parent_dir=parent_dir,
            original_dir=original_dir,
            ligand_id=ligand_id,
            env=env,
            cgenff_path=cgenff_path,
            protein_id=protein_id,
            rna=rna,
            raw_pdb_to_mol2=False,          # There is use for this as OpenBabel fails (silently) to treat certain structures e.g. triazoles, and creates a different molecule
            system_ph=7.4,                  # If a ligand is missing all hydrogens, it will be protonated by OpenBabel according to this pH
            ligand_input_sanitation=True,   # if False, ligands will not be checked for hydrogens, repeated atom names/numbers and double-uppercase element names
            relax_input_segid_filter=False, # if False, only PROA and HETA will be used
            #include_xray_water=False,      # Whether to include x-ray water molecules
            segid_filter=[                  # More precise control over segids included in waterboxes/complexes
            "PROB",
            "PROC",
            "PROD",
            "HETB",
            "HETC",
            "HETD",
            "WATA",
            "WATB",
            "WATC",
            "XRDA", # The XRDx segids are made by _add_segids()
            "XRDB", # The XRDx segids are made by _add_segids()
            "XRDC", # The XRDx segids are made by _add_segids()
            "XRDD", # The XRDx segids are made by _add_segids()
            "SOLV",
            "IONS",
        ],
        )
        
        # Create CHARMM Coordinate files
        used_segids = preparation.createCRDfiles()
        logging.info(
            f"Input parsing finished. Segids {', '.join(used_segids)} will be used "\
            f"for CRD generation and related tasks"
        )

        ########################################################################
        
        # Edit the CHARMM-GUI scripts
        # Instantiate class
        charmmManipulation = CharmmManipulation(
            parent_dir=parent_dir,
            ligand_id=ligand_id,
            original_dir=original_dir,
            het_resnames=preparation.het_resnames,
            env=env,
        )
        
        # Copy Files from the template folder
        charmmManipulation.copyFiles()

        # Modify step1_pdbreader.inp to read in correct amount of chains/residues
        charmmManipulation.modifyStep1(used_segids)

        # Run Charmm giving the correct executable path
        charmmManipulation.executeCHARMM(charmm_exe="charmm")
        
        # Convert the CHARMM system to an OpenMM system
        charmmManipulation.createOpenMMSystem()
        
        # Apply HMR to the OpenMM system
        charmmManipulation.applyHMR()
        
        # Create YAML file for ASFE simulations using TF (requires ./parent/config directory)
        charmmManipulation.createTFYamlFile(dt=0.002, nstep=2500000)
        
        # Space to pause and consider
        logging.info(f"{str(' '):>90}")
        
