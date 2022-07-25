#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:50:05 2022

@author: Johannes Karwanoupoulos and Åsmund Kaupang
"""
import sys
import os
import glob
import shutil
import subprocess

################################################################################
### FUNCTIONS
################################################################################
# HET BLOCK GENERATION
################################################################################
def genPROblock(protein_id, segment_id="PROA"):
    
    block = f"!-------------------------------------------------------------------------------\n"\
            f"! PROTEIN:\n"\
            f"open read card unit 10 name {protein_id.lower()}_{segment_id.lower()}.crd\n"\
            f"read sequence coor card unit 10 resid\n"\
            f"generate {segment_id.upper()} setup warn first NTER last CTER\n"\
            f"\n"\
            f"open read unit 10 card name {protein_id.lower()}_{segment_id.lower()}.crd\n"\
            f"read coor unit 10 card resid\n"\
            f"!-------------------------------------------------------------------------------\n"

    return block

def genHETblock(ligand_id, segment_id="HETA"):
    
    block = f"!-------------------------------------------------------------------------------\n"\
            f"open read card unit 10 name {ligand_id.lower()}.crd\n"\
            f"read sequence coor card unit 10 resid\n"\
            f"generate {segment_id.upper()} setup warn first none last none\n"\
            f"\n"\
            f"open read unit 10 card name {ligand_id.lower()}.crd\n"\
            f"read coor unit 10 card resid\n"\
            f"!-------------------------------------------------------------------------------\n"

    return block

def genWATblock(protein_id, segment_id="WATA"):
    
    block = f"!-------------------------------------------------------------------------------\n"\
            f"open read card unit 10 name {protein_id.lower()}_{segment_id.lower()}.crd\n"\
            f"read sequence coor card unit 10 resid\n"\
            f"generate {segment_id.upper()} setup warn noangle nodihedral\n"\
            f"\n"\
            f"open read unit 10 card name {protein_id.lower()}_{segment_id.lower()}.crd\n"\
            f"read coor unit 10 card resid\n"\
            f"!-------------------------------------------------------------------------------\n"

    return block

################################################################################
# HBUILD CONTROL
################################################################################
def hbuild_preserve_explicit_H(segids):

    segid_head = f"hbuild sele hydr .and. .not. ("
    segid_vars = ""
    for idx, segid in enumerate(segids):
        if idx == len(segids) - 1:
            segid_add = f"segid {segid}"
        else:
            segid_add = f"segid {segid} .or. "
        segid_vars = segid_vars + segid_add
    segid_tail = f") end"

    block = f"!-------------------------------------------------------------------------------\n"\
            f"! MOD: HBUILD control - preserve explicit H-coordinates\n"\
            f"prnlev 5\n"\
            f"echo START_HBUILD\n"\
            f"{segid_head + segid_vars + segid_tail}\n"\
            f"echo END_HBUILD\n"\
            f"!-------------------------------------------------------------------------------\n"
    return block

################################################################################
# PRINT USED PARAMETERS
################################################################################
def insert_parameter_print():
    block = f"!-------------------------------------------------------------------------------\n"\
            f"! MOD: Print used parameters\n"\
            f"echo START_PAR\n"\
            f"print para used\n"\
            f"echo END_PAR\n"\
            f"!-------------------------------------------------------------------------------\n"
    return block

################################################################################
# OTHER FUNCTIONS
################################################################################    

def streamWriter(work_dir, stream_name, blocks_as_lst):
    with open(f"{work_dir}/{stream_name}.str", 'w') as ostr:
        for block in blocks_as_lst:
            ostr.write(block)

def inputFileInserter(inpfile, cases, blocks, inversions):
    assert len(cases) == len(blocks)

    # Check for the presence of a backup, and start from this
    # or make a backup if none exists
    if os.path.exists(f"{inpfile}.premod"):
        shutil.copy(f"{inpfile}.premod", f"{inpfile}")
    else:
        shutil.copy(f"{inpfile}", f"{inpfile}.premod")

    inp_backup = f"{inpfile}.premod"
    outfile = inpfile

    with open(outfile, 'w') as out:
        with open(inp_backup, 'r') as inf:
            for line in inf:
                for case, block, inversion in zip(cases, blocks, inversions):
                    if case in line:
                        if inversion == True:
                            line = f"{block}\n"\
                                   f"\n"\
                                   f"{case}\n"
                        elif inversion == None:
                            line = f"{block}\n"
                        else:
                            line = f"{case}\n"\
                                   f"\n"\
                                   f"{block}\n"
                out.write(line)

def makeFolder(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
            
def makeTFFolderStructure(ligand_id, parent_dir="."):
    makeFolder(f"{parent_dir}/{ligand_id}")
    makeFolder(f"{parent_dir}/{ligand_id}/complex")
    makeFolder(f"{parent_dir}/{ligand_id}/waterbox")
    
def getTopparFromLocalCGenFF(ligands_dir, ligand_id, ligand_ext="mol2", cgenff_path=False, parent_dir="."):
    cgenff_bin = None
    cgenff_output = None

    
    # If no particular path is given, check whether CGenFF is available
    if cgenff_path == False:
        cgenff_path = shutil.which('cgenff')
        if cgenff_path == None:
            print("This function requires cgenff.")
            print("Please install it in the active environment or point the routine")
            print("to the right path using the key cgenff_path='/path/to/cgenff' .")
        else:
            cgenff_bin = cgenff_path
    else:
        cgenff_bin = cgenff_path

    # CGenFF exists - start program
    if cgenff_bin != None:
        
        # Run CGenFF
        ligand_path = f"{parent_dir}/{ligands_dir}/{ligand_id}.{ligand_ext}"
       
        cgenff_output = subprocess.run(#
                                        [cgenff_bin] + [ligand_path] + ["-v"] + ["-f"] + [f"{ligand_id}/{ligand_id}.str"] + ["-m"] + [f"{ligand_id}/{ligand_id}.log"],#
                                        #[cgenff_bin] + [ligand_path] + ["-v"] + ["-m"] + [f"{ligand_id}.log"],#
                                        capture_output=True,
                                        text=True#
                                        )
        
        # Evaluate the subprocess return code
        if cgenff_output.returncode == 1:
            print(f"CGenFF returned an error after being called with:\n")
            print(' '.join(cgenff_output.args))
            print(cgenff_output.stdout)
            print(cgenff_output.stderr)
        else:
            print(f"CGenFF executed successfully")
            #print(cgenff_output.stdout)
            #print(cgenff_output.stderr)
            
    return cgenff_output

