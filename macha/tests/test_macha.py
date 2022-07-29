
import pytest
from ..functions import check_ligands, Preparation, CharmmManipulation



def test_import():
    import parmed as pm


def test_preparation():

    ligand_id = "../data/templates/cdk2_32.pdb"
    parent_dir = "."
    original_dir = "../data/original"

    preparation = Preparation(parent_dir=parent_dir, ligand_id=ligand_id,original_dir=original_dir,env="waterbox")

    # Make a Transformato style folder structure below a folder bearing
    # the name of the ligand
    preparation.makeTFFolderStructure()
    segids = preparation.createCRDfiles()