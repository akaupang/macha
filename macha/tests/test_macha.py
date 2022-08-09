from random import randint
import pytest
from macha.functions import Preparation, checkInput
import os
import parmed as pm

ligand_id = "../data/original/cdk2_32.pdb"
parent_dir = "."
original_dir = "macha/data/original"

def test_import():
    import parmed as pm
    import openbabel

def test_checkInput():

    proteinID,ligandID = checkInput(parent_dir=parent_dir, original_dir=original_dir, protein_name=None, input_ext="pdb")

    # assure that all ligands are found by the function
    if not proteinID:
        assert len(os.listdir(original_dir)) == len(ligandID)

def test_createFolders():

    preparation = Preparation(parent_dir=parent_dir, ligand_id="cdk2_l32",original_dir=original_dir,env="waterbox")
    pdb_macha = preparation.pdb_file
    pdb_orig = pm.load_file("macha/data/original/cdk2_l32.pdb")
   
    number = randint(0,9065)
    assert str(pdb_macha.atoms[number]) == str(pdb_orig.atoms[number])

