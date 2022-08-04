import pytest
from macha.functions import Preparation, check_ligands


ligand_id = "../data/templates/cdk2_32.pdb"
parent_dir = "."
original_dir = "../data/original"

def test_import():
    import parmed as pm
    import openbabel


def test_createFolders():

    preparation = Preparation(parent_dir=parent_dir, ligand_id=ligand_id,original_dir=original_dir,env="waterbox")
    preparation.makeTFFolderStructure()