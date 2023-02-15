#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:50:05 2022

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""

import sys
import os
import glob
import shutil
import subprocess
import parmed as pm
from openbabel import openbabel
import string
import numpy as np
from .charmm_factory import CharmmFactory

import warnings

# Supress ParmEd warnings
warnings.filterwarnings("ignore", module="parmed")


def checkInput(
    parent_dir=".", original_dir="original", protein_name=None, input_ext="pdb"
):

    # Make sure the input directory exists
    if not os.path.isdir(f"{parent_dir}/{original_dir}"):
        sys.exit(f"Input directory: {parent_dir}/{original_dir} not found!")

    # Look for an input protein for ligand exchange
    if protein_name != None:
        protein_path = f"{parent_dir}/{original_dir}/{protein_name}.{input_ext}"
        if not os.path.isfile(protein_path):
            print(f"No protein for ligand exchange found.")
            protein_id = None
        else:
            print(f"Protein for ligand exchange found at: {protein_path}")
            protein_id = os.path.splitext(os.path.basename(protein_path))[0]
            print(f"Protein ID: {protein_id}")
    else:
        protein_id = None

    if protein_id == None:
        print(
            """
        No protein was specified for ligand exchange with multiple ligands.
        The input is thus assumed to consist either of complexes from which to create
        waterboxes + complexes, OR of single ligands from which to create waterboxes
        (if the latter is the case, use the -nc switch for efficiency).
        """
        )

    # Create the master list of the ligand/system names and fill it
    ligand_ids = []
    # Add multiple ligands
    for ifile in glob.glob(f"{parent_dir}/{original_dir}/*.{input_ext}"):
        lig_id = os.path.splitext(os.path.basename(ifile))[0]
        # But don't add the protein/complex
        if lig_id == protein_id:
            pass
        else:
            ligand_ids.append(lig_id)  #
    return protein_id, ligand_ids

###############################################################################
class Preparation:
###############################################################################
    def __init__(
        self,
        parent_dir,
        original_dir,
        ligand_id,
        env,
        protein_id=None,
        small_molecule=False,
        rna=False,
        system_ph=7.4,
        input_sanitation=True,
    ):

        """
        This class prepares everything for further use with CHARMM. The PDB
        files are sliced into pieces and the ligand is converted to a mol2 file.
        A local version of CGenFF creates a stream file for the ligand.
        """
        self.parent_dir = parent_dir
        self.original_dir = original_dir
        self.ligand_id = ligand_id
        self.protein_id = protein_id
        self.resname = str
        self.env: str = env
        self.small_molecule: bool = small_molecule
        self.rna: bool = rna
        self.input_type = None
        self.system_ph = system_ph        
        
        # Distinguished treatment if the current ligand is to be merged with a protein
        if self.protein_id == None:
            # Original input file
            self.orig_ligand_input = (
                f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb"
            )
            
            # Optionally disabled input sanitation
            if input_sanitation == True:
                self.sanitizeInput(self.orig_ligand_input)
            
            # ENVIRONMENTS
            # For single ligands -> waterbox
            # Perhaps also useful for:
            # complex.pdb -> complex + ligand (from complex) in waterbox
            # Load the PDB file into ParmEd
            self.pm_obj = pm.load_file(self.orig_ligand_input, structure=True)
            
            if self.env == "single_strand":
                self.pm_obj = pm.load_file(self.orig_ligand_input, structure=True)
                self.pm_obj = self.pm_obj["A", :, :]    # select only CHAIN A 
                
        else:
            # Original input file
            self.orig_ligand_input = (
                f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb"
            )
            
            # Optionally disabled input sanitation
            if input_sanitation == True:
                self.sanitizeInput(self.orig_ligand_input)

            self.orig_protein_input = (
                f"{self.parent_dir}/{self.original_dir}/{self.protein_id}.pdb"
            )
            # Input proteins are not sanitized!
                        
            # ENVIRONMENTS            
            if self.env == "waterbox" or self.env == "double_strand":
                # Run like normal single ligand
                self.pm_obj = pm.load_file(self.orig_ligand_input, structure=True)
                
            elif self.env == "single_strand":
                self.pm_obj = pm.load_file(self.orig_ligand_input, structure=True)
                self.pm_obj = self.pm_obj["A", :, :]    # select only CHAIN A

            elif self.env == "complex":
                # Start the merger function
                self.mergeToComplex() 

            else:
                sys.exit(f"Unrecognized environment: {self.env}")
    
    def makeTFFolderStructure(self):
        self._make_folder(f"{self.parent_dir}/config")
        self._make_folder(f"{self.parent_dir}/{self.ligand_id}")
        self._make_folder(f"{self.parent_dir}/{self.ligand_id}/{self.env}")
        self._make_folder(f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm")
        self._make_folder(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/restraints"
        )
        
    def _make_folder(self, path):
        try:
            os.makedirs(path)
            print(f"Creating folder {path} for the {self.env}")
        except OSError:
            print(f"The folder {path} for the {self.env} exists - we will use it")
            if not os.path.isdir(path):
                raise
            
    def sanitizeInput(self, structure_file_path):
        import re
        
        # This only works for PDBs
        assert structure_file_path.endswith(".pdb")
        
        # Backup original file if a backup does not exist
        if os.path.isfile(f"{structure_file_path}.bck"):
            shutil.copy(
                f"{structure_file_path}.bck",
                structure_file_path,
            )
        else:
            shutil.copy(
                structure_file_path,
                f"{structure_file_path}.bck",
            )
        
        # Compile regular expressions
        UU_Case = re.compile('([A-Z][A-Z])')
                       
        # Define private functions
        def _convert_double_uppercase_to_upperlower(match_obj):
            if match_obj.group(1) is not None:
                return match_obj.group(1)[0].upper() + match_obj.group(1)[1].lower()
        
        
        # Open input and output, and write eventual changes
        # This checks and fixes:
        # - double uppercase atom names
        # - instances of duplicate atom name+number
        #
        with open(f"{structure_file_path}.bck", "r") as isf:
            with open(structure_file_path, "w") as osf:
                atom_names = []
                for line in isf:
                    if ((line.startswith("HETATM")) or (line.startswith("ATOM"))):
                        
                        # Gather portions of the PDB file
                        linestart=line[0:12]
                        atomname=line[12:16]
                        linemiddle=line[16:76]
                        elementname=line[76:]
                             
                        # Check and fix double uppercase
                        if re.match(UU_Case, atomname):
                            UU_an = atomname.replace(' ','')
                            atomname = re.sub(UU_Case, _convert_double_uppercase_to_upperlower, atomname)
                            Uu_an = atomname.replace(' ','')
                            print(
                                f"Double uppercase atom name {UU_an} "\
                                f"replaced with {Uu_an}"
                            )
                        
                        if re.match(UU_Case, elementname):
                            UU_en = elementname.replace(' ','').strip('\n')
                            elementname = re.sub(UU_Case, _convert_double_uppercase_to_upperlower, elementname)
                            Uu_en = elementname.replace(' ','').strip('\n')
                            print(
                                f"Double uppercase atom name {UU_en} "\
                                f"replaced with {Uu_en}"
                            )
                            
                        # Check if atomname has been seen before
                        this_an = atomname.replace(' ','')
                        if this_an in atom_names:
                            # Get the existing atom name - there should be only one
                            exist_an = [name for name in atom_names if name == this_an][0]
                            
                            # Split the atom name and the atom number
                            nn_lst = re.split('(\d+)', exist_an)
                            
                            # Add 1 to the atom number and recompile atom name/number
                            new_aname = f"{nn_lst[0]}"
                            new_anumber = f"{int(nn_lst[1])+1}"
                            
                            # Format the new atomname portion for the line
                            atomname = f"{new_aname:>2}{new_anumber:<2}" # 12:16 = 4 characters
                                                        
                            # Add modified atomname to the list, in case duplication propagates
                            atom_names.append(new_aname + new_anumber)
                            
                            # Let the user know that this replacement happened
                            print(f"Atom name {exist_an} exists. Replaced with {new_aname + new_anumber}")

                        else:
                            # Add unique atomname to a list for comparison
                            atom_names.append(this_an)

                        # Write the finished line
                        osf.write(linestart + atomname + linemiddle + elementname)
                                                    
                    else:
                        osf.write(line)

        




        
        


    def mergeToComplex(self):
        print("* * * Begin merge to complex * * *")
        # Here, it is practical to borrow the "self.pm_obj" temporarily
        # for both the input protein and later for the ligand. This solution
        # benefits from the fact that most of the internal functions in this
        # class already work on this object. The merge function will in the end
        # restore "self.pm_obj" as the complex.

        # Load the protein
        print(f"Processing protein (input typing, segid addition, HETx removal)")
        self.pm_obj = pm.load_file(
            self.orig_protein_input, structure=True
        )

        # What kind of protein PDB file has been provided?
        # Remember: segids will be added if there are none
        segids, pm_obj_df = self.checkInputType(is_complex=True)

        # Get the parts of the ParmEd object (here, the protein) that does not 
        # contain "HET" in the segid column. No segids other than HET are excluded!
        # Note the tilde to invert the boolean mask!
        # Copy the resulting ParmEd object
        apo_protein_pm_obj = self.pm_obj[~pm_obj_df.segid.str.contains("HET")]

        # Load the ligand, overwriting the common ParmEd object
        # HERE THE LIGAND WITH HYDROGENS FROM OPENBABEL MUST BE LOADED
        print(f"Processing ligand: {self.ligand_id} (input typing, segid addition)")
        self.pm_obj = pm.load_file(
            self.orig_ligand_input, structure=True
        )
        
        # Add segids to ligand ParmEd object if there are none
        segids, pm_obj_df = self.checkInputType(is_complex=False)

        # Copy the ligand ParmEd object
        ligand_pm_obj = self.pm_obj

        # Overwrite the self.pm_obj by joining the copied protein and ligand
        # ParmEd objects.
        # NOTE: No reordering to a typical PROA, HETA, IONS, WAT/SOLV order
        # is performed
        self.pm_obj = apo_protein_pm_obj + ligand_pm_obj

        print("* * * End merge * * *")

    def checkInputType(self, is_complex=False):
        # Since this function runs right after the __init__ function,
        # ParmEd objects of several origins need to be handled:
        # Ligands (or RNA) for waterboxes and ligands during mergetocomplex
        # Proteins during mergetocomplex
        # Complexes after mergetocomplex (actually redundant, since both ran through before the merge)
        
        # We will cater for two types of PDB files
        # 1: Canonical/MAESTRO-derived PDB files
        # 2: CHARMM-derived PDB files
        #
        # This distinction is made by looking for segids.
        #
        # PDB TYPE                  CHAIN IDs       SEGIDs
        # Canonical/MAESTRO         yes             no (element names)
        # CHARMM PDB                no              yes
        # CHARMM PDB ("official")   yes             yes
        
        # TREATMENT FOR ALL OBJECTS
        # Remove lone pairs (if any)
        self.pm_obj = self._remove_lp()

        # Check if the user has set the protonation states of ionizable residues
        # using either CHARMM or AMBER residue names (not quite comprehensively defined)
        # If not, ATTENTION, all occurrences of HIS will be renamed to HSD
        self.pm_obj = self._check_ionizable()

        # Store the ParmEd object as a dataframe for later use
        pm_obj_df = self.pm_obj.to_dataframe()        
        
        # THIS MAY BE THE PLACE TO DO BASIC INPUT SANITATION SUCH AS CHECKS FOR:
        # - doubly capitalized atom names, e.g. CL, NA
        # - atoms with the same name

        # # Count the segids in the ParmEd object (if any)
        segids = set(i.residue.segid for i in self.pm_obj)

        # For RNA (some extra treatment is necessary)
        if self.rna:
            self.pm_obj = self._add_segids_rna(pm_obj_df)
            self.pm_obj = self._create_tlc_rna()
            segids = set(i.residue.segid for i in self.pm_obj)
            pm_obj_df = self.pm_obj.to_dataframe()
        
        # For canonical/MAESTRO-derived PDB files containing chain IDs
        elif segids == {""}:  # empty set of segids
            print(f"Processing a ParmEd object without segids (e.g. from a canonical/Maestro based PDB file)")
            # Add segids to the Parmed object (NEEDED FOR CHECK HYDROGEN FUNCTION)
            self.pm_obj = self._add_segids(pm_obj_df)
            
            # Update the dataframe and segid list
            pm_obj_df = self.pm_obj.to_dataframe()     
            segids = set(i.residue.segid for i in self.pm_obj)
            
            # Check for hydrogens and update the dataframe
            self._check_for_hydrogens(is_complex=is_complex)
            pm_obj_df = self.pm_obj.to_dataframe()     

            # Add segids to the updated Parmed object
            self.pm_obj = self._add_segids(pm_obj_df)
            
            # Update the dataframe before it is returned
            pm_obj_df = self.pm_obj.to_dataframe()
            
        # For CHARMM(-GUI)-generated PDB files
        else:
            print(f"Processing a ParmEd object with segids (e.g. from a CHARMM PDB file)")
            # Check for hydrogens and update the dataframe
            self._check_for_hydrogens(is_complex=is_complex)
            pm_obj_df = self.pm_obj.to_dataframe()     
            
            # Add segids to the Parmed object in case it was changed
            self.pm_obj = self._add_segids(pm_obj_df)
            
            # Update the dataframe before it is returned in case it was changed
            pm_obj_df = self.pm_obj.to_dataframe()            
        
        return segids, pm_obj_df

    def _remove_lp(self):
        # Will check if there are lone pairs and remove them
        # CGenFF will add them later on
        lps = []
        for atom in self.pm_obj:
            if atom.name.upper().startswith("LP"): # Check for lp, Lp, lP, LP
                print(f"We will remove lonepair {atom.idx} {atom.name}, {atom}")
                lps.append(atom.idx)

        for i in range(len(lps)):
            self.pm_obj.strip(f"@{lps[i]+1-i}")

        return self.pm_obj

    def _check_ionizable(self):
        # A preliminary list of residue names for titrable residues that indicate 
        # that the user has set them manually appears below
        # AMBER names will be converted to CHARMM names, no questions asked.
        # TODO This list should be expanded to be comprehensive
        ionized_res = [
            "HID",
            "HIE",
            "HIP",
            "HSD",
            "HSE",
            "HSP",
            "CYX",
            "CYM",
            "ASPP",
            "GLUP",
        ]

        for (
            res
        ) in (
            self.pm_obj.residues
        ):  
            if res.name in ionized_res:
                # Harmonize to CHARMM names (from AMBER names)
                # TODO This section should be expanded to be comprehensive
                if res.name == "HID":
                    res.name = "HSD"
                elif res.name == "HIE":
                    res.name = "HSE"
                elif res.name == "HIP":
                    res.name = "HSP"
                elif res.name == "CYX":
                    res.name = "CYM"
                else:
                    pass

            elif res.name == "HIS":
                # TODO A lookup using PROPKA or another solution could be of interest here
                res.name = "HSD"  # ATTENTION!! Here we make all HIS to HSD

        return self.pm_obj

    def _check_for_hydrogens(self, is_complex=False):
        # Get the segids and dataframe
        segids = set(i.residue.segid for i in self.pm_obj)
        pm_obj_df = self.pm_obj.to_dataframe()
        
        # Normal run for ligands, etc.
        if is_complex == False:
            # Exclude segments that do not need explicit hydrogens (that ICBuild/HBuild can handle)
            het_segids = [segid for segid in segids if segid[:3] not in ["PRO", "WAT", "ION", "XRD"]]
            
            # THIS IS A LITTLE CLUNKY AT THE MOMENT, BUT MAY BE ADAPTED LATER FOR MULTIPLE LIGANDS*
            
            
            if len(het_segids) == 1:
                #for segid in het_segids:*
                segid = het_segids[0]
                struct = self.pm_obj[pm_obj_df.segid == segid]
                
                # Make sure there is only one residue in this segid*
                assert len(struct.residues) == 1
                
                # Work on the single residue in this segid*
                res = struct.residues[0]
                
                # Check for the presence of hydrogens in the structure                 
                if len([atm for atm in res.atoms if atm.element_name == 'H']) == 0:
                    print(f"No hydrogens found in residue {res.name} from segid {segid}")
                    
                    # Slice a ParmEd object based on the residue (in the segid) 
                    # that does not contain any hydrogens
                    res_obj = self.pm_obj[(pm_obj_df.segid == segid) | (pm_obj_df.resname == res.name)]
                    
                    # Set self.resname
                    self.resname = res.name
                    
                    # Make a folder for the ligand toppar files
                    self._make_folder(f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}")

                    # Save the object to a PDB file
                    res_obj.save(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_noH_PE.pdb",
                        overwrite=True,
                        )
                    
                    # Define an OpenBabel object
                    obConversion = openbabel.OBConversion()
                    mol = openbabel.OBMol()
                                    
                    # Read a ligand PDB file to OpenBabel object mol
                    obConversion.ReadFile(
                            mol,
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_noH_PE.pdb",
                        )   
                            
                    # ATTENTION: OpenBabel may or may not get this right
                    # May help if there is trouble:
                    # mol.UnsetFlag(openbabel.OB_PH_CORRECTED_MOL)
                    # mol.SetAutomaticFormalCharge(True)
                    mol.CorrectForPH(self.system_ph)
                    mol.AddHydrogens()
                    
                    # Save the hydrogenated OpenBabel object to a PDB file (unfortunately without resname and chain ID)        
                    obConversion.WriteFile(
                        mol,
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_withH_OB.pdb",
                    )
                    
                    # TAKE THE FIRST HALF OF THE NEW FILE FROM THE OLD (THE HEAVY ATOMS)      
                    with open(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_noH_PE.pdb",
                        "r"
                    ) as ipdb_noh:                         
                        with open(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_withH_FIX.pdb",
                            "w"
                        ) as opdb:
                            for line in ipdb_noh:
                                if line.startswith("END"):
                                    pass
                                else:
                                    opdb.write(line)  
                    
                    # THEN, TAKE THE HETATM/ATOM H LINES AND APPEND THESE AFTER THE HEAVY ATOMS, UPDATING MISSING INFORMATION
                    # ON THE FLY                                       
                    with open(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_withH_FIX.pdb",
                        "a"
                    ) as opdb:       
                        with open(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_withH_OB.pdb",
                            "r"
                        ) as ipdb_h:
                            h_num = 1
                            for line in ipdb_h:
                                if (((line.startswith("HETATM")) or (line.startswith("ATOM"))) and (line[13:16].strip(" ").startswith("H"))):
                                    line = line.replace("UNL", res.name)
                                    line = line.replace(f" H   {res.name}",f" H{h_num:<3}{res.name}")
                                    line = line.replace(f"{res.name}  ", f"{res.name} {res.chain}")
                                    opdb.write(line)
                                    h_num += 1
                                    
                            opdb.write(f"END\n")

                    # Remake the ParmEd object
                    self.pm_obj = pm.load_file(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_withH_FIX.pdb",
                        structure=True
                    )
                    
                    print(
                        f"Added "\
                        f"{len([atm for atm in self.pm_obj.atoms if atm.element_name == 'H'])}"\
                        f" hydrogens to residue {res.name} for pH {self.system_ph:.1f}"
                    )
                    
                    # Resave the ParmEd object to view how it was perceived
                    self.pm_obj.save(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_withH_FIX_PE.pdb",
                        overwrite=True,
                    )

                else:
                    print(f"Hydrogens found in residue {res.name} from segid {segid}")
            
            else:
                sys.exit("More than 1 HET segid found - multiple ligands are not supported!")
            
        # Or skip checks, in case of protein from complex
        else:
            print(
                f"Skip checking whether hydrogens are defined in a protein-ligand complex,\n"\
                f"since the apo protein will be produced in the next step (and all other segids are discarded)."
            )    
      
    




    def _add_segids_rna(self, pm_obj_df):

        chids = set(i.residue.chain for i in self.pm_obj)
        print(f"Found chain IDs: " + ", ".join(chids))
        for chain in chids:
            for res in self.pm_obj.view[
                pm_obj_df.chain == f"{chain}"
            ].residues:  # produce a select view of this chain's residues using a boolean mask
                res.segid = f"RNA{chain}"

        return self.pm_obj
    
    def _create_tlc_rna(self):

        for atom in self.pm_obj.atoms:
            if atom.residue.name == "G":
                atom.residue.name = "GUA"
            elif atom.residue.name == "C":
                atom.residue.name = "CYT"
            elif atom.residue.name == "A":
                atom.residue.name = "ADE"
            elif atom.residue.name == "U":
                atom.residue.name = "URA"
            elif atom.residue.name == "T":
                atom.residue.name = "THY"
            elif atom.residue.name == "I":
                atom.residue.name = "INO"

        return self.pm_obj

    def _add_segids(self, pm_obj_df):
        # This function adds segids to the Parmed object from a PDB that does
        # not contain segids (only chain ids), ensuring that segids can be
        # used in making of CRD files later.
        #
        # The function evaluates the residue names and gives those that are 
        # "normal" amino acids the segid PROx, where x is their chain ID.
        # Waters (HOH) are assigned the segid WAT+chain ID
        # Residue names in the "HET-exclusion" list below (e.g. common
        # additives in the crystallization buffer, ions, etc) are given
        # the segid XRD+chain ID.
        # Residues that do not correspond to normal amino acid names or
        # the names given in the exclusion list are given the segid HETx, 
        # in which x is a letter in ABCDEFGH...

        # TODO This list should be expanded to be comprehensive (look in CHARMM toppar for RESI and PRES)
        aa_res = [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TPO",
            "TRP",
            "TYR",
            "VAL",
            "HID",
            "HIE",
            "HIP",
            "HSD",
            "HSE",
            "HSP",
            "CYX",
            "CYM",
            "ASPP",
            "GLUP",
        ]  # we need to find the residue of the ligand which should be the only one being not an aa

        exclude_res = [# TODO Many more common crystallization additives should be added to this list
            "1PE"
        ]

        # ions_res =[] # TODO A list of ions could be added for their particular handling

        # RUDIMENTARY TEST FOR EXISTING SEGIDS (NEGATIVE)
        # TRUE IF NOT ALL ENTRIES HAVE A SEGID
        if len(pm_obj_df) != len([line for line in pm_obj_df.segid if line != ""]):
            # Make a list of ABCD... for use as x in HETx
            het_letters = list(string.ascii_uppercase)

            # TODO Make this loop handle multiple ligands in one chain
            chids = set(i.residue.chain for i in self.pm_obj)
            # DEBUG
            print(f"Found chain IDs: " + ", ".join(chids))
            
            for (
                chain
            ) in chids:
                # DEBUG
                #print(f"Working on chain {chain}")
        
                # Set a temporary previous residue number
                prev_resnum = -1

                # Loop over the residues in the current chain
                for res in self.pm_obj.view[
                    pm_obj_df.chain == f"{chain}"
                ].residues:  # produce a select view of this chain's residues using a boolean mask
                    # DEBUG
                    #print(f"Working on residue {res.name}")
                    resnum = res.number
                    if res.name in aa_res:
                        res.segid = f"PRO{chain}"
                    elif ((res.name == "HOH") or (res.name.startswith("TIP"))):
                        res.segid = f"WAT{chain}"
                    elif res.name in exclude_res:
                        res.segid = f"XRD{chain}"
                    else:
                        if resnum == prev_resnum:
                            pass
                        else:
                            res.segid = f"HET{het_letters.pop(0)}"
                    prev_resnum = resnum
        else:
            print(f"All entries appear to have a segid already")

        return self.pm_obj
    
    def createCRDfiles(self, segids, pm_obj_df):

        exclude_segids = [
            "HETB",
            "HETC",
            "SOLV",
            "IONS",
            "WATA",
            "WATB",
            "WATC",
            "XRDA",
            "XRDB",
            "XRDC",
        ]
        
        used_segids = []
        for segid in segids:
            if segid not in exclude_segids:  # multiple ligands are excluded for now
                # Store in the segids to be given to CharmmManipulation
                used_segids.append(segid)

                # WATERBOX ENVIRONMENT
                if self.env == "waterbox": #self.env != "complex": # 2023: why not == "waterbox" ?
                    # For ligands
                    if segid.startswith("HET"):
                        # DEBUG
                        #print(f"Env: {self.env}, Segid: {segid}")
                        
                        # Note the residue name for checks
                        self.resname = (
                            self.pm_obj[pm_obj_df.segid == f"{segid}"].residues[0].name
                        )
                        # resname should be a 3 or 4 letter code
                        assert len(self.resname) < 5

                        # Make a folder for the ligand files in the waterbox directory
                        self._make_folder(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}"
                        )
                        # Write the CHARMM CRD file from the ParmEd object
                        self.pm_obj[pm_obj_df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        
                        # If the switch is active, copy the original PDB file to the ligand folder
                        # This seems RISKY, since this PDB file has not been vetted by checkInputType 
                        # (e.g for hydrogen existence)
                        if self.small_molecule:
                            # We copy the pdf file since there shouldn't be any changes
                            shutil.copy(
                                f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb",
                                f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                            )
                        else:
                            # If the switch is not active, 
                            # save the CHARMM PDB file of the ligand from the ParmEd object
 
                            self.pm_obj[pm_obj_df.segid == f"{segid}"].save(
                                f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                                overwrite=True,
                            )
                    
                    # For RNA, a segid-based naming applies (instead of resname-based as above)
                    elif segid.startswith("RNA"):
                        # DEBUG
                        #print(f"env: {self.env}, segid: {segid}")
    
                        self.pm_obj[pm_obj_df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )

                # COMPLEX ENVIRONMENT
                elif self.env == "complex":
                    if segid.startswith("HET"):
                        # DEBUG
                        #print(f"Env: {self.env}, Segid: {segid}")

                        # TODO: Should this check be more generally applied?
                        # Note the residue name for checks
                        self.resname = (
                            self.pm_obj[pm_obj_df.segid == f"{segid}"].residues[0].name
                        )
                        # resname should be a 3 or 4 letters code
                        assert len(self.resname) < 5

                        # Make a folder for the ligand in the complex directory
                        self._make_folder(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}"
                        )
                        
                        # Save the CHARMM CRD file of the ligand from the ParmEd object
                        self.pm_obj[pm_obj_df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        
                        # Save the CHARMM PDB file of the ligand from the ParmEd object
                        self.pm_obj[pm_obj_df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                            overwrite=True,
                        )

                    else: 
                        # DEBUG
                        #print(f"Env: {self.env}, Segid: {segid}")

                        # Save the protein, as well as eventual other segids,
                        # as CHARMM CRD and -PDB files, directly in the complex directory
                        self.pm_obj[pm_obj_df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                            overwrite=True,
                        )
                        self.pm_obj[pm_obj_df.segid == f"{segid}"].save(
                            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.pdb",
                            overwrite=True,
                        )
                else:
                    # Alert if unsupported environments are requested or if input has typos
                    sys.exit(f"Unrecognized environment: {self.env}")

        return segids, used_segids
    
    def getTopparFromLocalCGenFF(
        self,
        cgenff_path=False,
    ):
        cgenff_bin = None
        cgenff_output = None

        # If no particular path is given, check whether CGenFF is available
        if cgenff_path == False:
            cgenff_path = shutil.which("cgenff")
            if cgenff_path == None:
                print(f"This function requires cgenff.\n"\
                      f"Please install it in the active environment or point the routine\n"\
                      f"to the right path using the key cgenff_path='/path/to/cgenff' ."
                )
            else:
                cgenff_bin = cgenff_path
        else:
            cgenff_bin = cgenff_path

        # CGenFF exists - start program
        if cgenff_bin != None:
            
            # CGenFF needs a mol2 file as input file
            self._create_mol2_sdf_file()
            
            # If a ligand toppar stream file exists (from a previous rung), 
            # we will remove it to prevent CGenFF from appending the new topology
            # and parameters at the end of the old one
            stream_file = f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.str"
            if os.path.isfile(stream_file):
                os.remove(stream_file)

            # Run CGenFF
            ligand_path = f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}"

            cgenff_output = subprocess.run(
                [cgenff_bin]
                + [f"{ligand_path}.mol2"]
                + ["-v"]
                + ["-f"]
                + [f"{stream_file}"]
                + ["-m"]
                + [f"{ligand_path}.log"],
                capture_output=True,
                text=True,
            )

            # Evaluate the subprocess return code
            if cgenff_output.returncode == 1:
                print(f"CGenFF returned an error after being called with:\n"\
                      f"{' '.join(cgenff_output.args)}\n"\
                      f"{cgenff_output.stdout}\n"\
                      f"{cgenff_output.stderr}\n"
                )
            else:
                print(f"CGenFF executed successfully")

            # Modify the resname in the new stream file
            self._modify_resname_in_stream()

        return cgenff_output
    
    def _create_mol2_sdf_file(self):

        print(f"Converting ligand {self.ligand_id} to MOL2 and SDF files")
        obConversion = openbabel.OBConversion()
        mol = openbabel.OBMol()
                
        # Read a ligand PDB file to OpenBabel object mol
        
        # THIS APPEARS TO BE REDUNDANT SINCE small_molecule == True already copies 
        # the original PDB to the ligand folder
        
        # # From the "original" directory
        # if self.small_molecule:
        #     obConversion.ReadFile(
        #         mol,
        #         f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb",
        #     )
        # # From the "ligandID/environment/residueName" directory, e.g. "cpd10/waterbox/dy6/dy6.pdb"
        # else:
        #     obConversion.ReadFile(
        #         mol,
        #         f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
        #     )
        
        obConversion.ReadFile(
                mol,
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
            )
        
        assert (mol.NumResidues()) == 1
        
        # Write new files
        # Convert mol from PDB to MOL2
        obConversion.SetInAndOutFormats("pdb", "mol2")
        obConversion.WriteFile(
            mol,
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.mol2",
        )
        
        # Convert mol from PDB to SDF
        obConversion.SetInAndOutFormats("pdb", "sdf")
        obConversion.WriteFile(
            mol,
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.sdf",
        )
        
    def _modify_resname_in_stream(self):

        # Open the ligand topology and parameter stream file produced by CGenFF
        stream_in = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.str",
            "rt",
        )
        # Open a new ligand topology and parameter stream file
        stream_out = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}_tmp.str",
            "wt",
        )
        
        # Correct the residue name in the stream file
        for line in stream_in:
            if line.startswith("RESI"):
                stream_out.write(line.replace(line.split()[1], self.resname))
            else:
                stream_out.write(line)

        stream_in.close()
        stream_out.close()

        shutil.copy(stream_out.name, stream_in.name)

    
###############################################################################
class CharmmManipulation:
###############################################################################
    def __init__(
        self,
        parent_dir,
        ligand_id,
        original_dir,
        resname,
        env,
        include_ions=True,
        ion_name="POT", # EXPAND FOR CATION AND ANION?
        ion_conc=0.15,
        default_path="",
    ):
        """
        CHARMM related files like the toppar file are modified, later the CHARMM executable is
        executed
        """
        self.parent_dir: str = parent_dir
        self.ligand_id: str = ligand_id
        self.original_dir: str = original_dir
        self.resname: str = resname
        self.env: str = env
        self.include_ions: bool = include_ions
        self.ion_name: str = ion_name
        self.ion_conc: float = ion_conc
        if not default_path:
            self.default_path: str = self.get_default_path()
        else:
            self.default_path = default_path

    def get_default_path(self):

        import macha

        return f"{macha.__path__[0]}/data/templates/default/"
    
    def copyFiles(self):
        # Copy files from template path
        # CHARMM related files
        for file in glob.glob(f"{self.default_path}/*[!omm_*][!toppar][!__pycache__]*"):
            shutil.copy(file, f"{self.parent_dir}/{self.ligand_id}/{self.env}/")

        # CHECK FOR USER INPUT ON ION TYPE AND ION STRENGTH
        # If someething other than the default is requested, edit the step2.2 ions 
        # stream file on the fly.
        if self.ion_name != "POT" or self.ion_conc != 0.15:
            self._modifyIonsCountFile()

        # CheckFFT script
        shutil.copy(
            f"{self.default_path}/checkfft.py",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/checkfft.py",
        )

        # The CHARMM toppar.str (which later will be converted to OpenMM format)
        shutil.copy(
            f"{self.default_path}/toppar.str",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar.str",
        )
        # The whole toppar directory
        try:
            shutil.copytree(
                f"{self.default_path}/toppar",
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar",
            )
        except:
            print(f"Toppar directory is already available")
    
    def _modifyIonsCountFile(self):

        # Create a target for a new step2.2 ions stream file
        fout = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step2.2_ions_count_tmp.str",
            "wt",
        )
        
        # Read the existing template step2.2 ions stream file and replace the occurrences
        # of "POT" and "0.15", with the cation type and concentration requested by the user
        with open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step2.2_ions_count.str",
            "r+",
        ) as f:
            for line in f:
                if "set pos" in line:
                    new_line = line.replace("POT", str(self.ion_name))
                    fout.write(new_line)
                elif "set conc" in line:
                    new_line = line.replace("0.15", str(self.ion_conc))
                    fout.write(new_line)
                else:
                    fout.write(line)
        fout.close()

        shutil.copy(
            fout.name,
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step2.2_ions_count.str",
        )
        os.remove(fout.name)

        pass
    
    def modifyStep1(self, segids):

        # NOTE The function _appendLigandTopparToTopparStream (3 lines) is only used here
        # Perhaps it is better to include it directly?
        # Also, it is unclear what throws the TypeError exception below
        # - this should be commented on to explain the integration of the
        # RNA "layer"
        
        # Append the ligand toppar to the CHARMM toppar stream
        try:
            self._appendLigandTopparToTopparStream()
        except TypeError:
            print("Won't add a ligand stream file - assuming that the ligand is an RNA strand")


        # Open a a target for the new input file
        fout = open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step1_pdbreader_tmp.inp",
            "wt",
        )
        # Read the template input file and 

        header_section_active = False # A switch to mark whether we are in the header 
                                      # section of step1_pdbreader as we read the file
                                      # and replace/copy its lines
        with open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step1_pdbreader.inp", "r+"
        ) as f:
            for line in f:
                if line.startswith("! Read PROA"): # Marks the start of the header section
                                                   # with the read/generate statements
                    header_section_active = True
                    f.readline() # What does this do?
                    
                    # Create read/generate statements from segids
                    # NOTE WHY IS THIS NOT A FUNCTION IN THE CHARMMMANIPULATION CLASS, 
                    # BUT RATHER IN ITS OWN CLASS (CHARMMFACTORY)?
                    header_block = CharmmFactory.createHeader(segids, self.env)
                    fout.write(f"{header_block}\n")                    
                    
                if line.startswith("!Print heavy atoms with "): # Marks the end of the header section
                    header_section_active = False

                # If we are in the header section, we should not record the lines in the
                # template file since these have been replaced by new ones
                if header_section_active == True:
                    pass
                else: # The rest of the step1_pdbreader.inp however, we copy verbatim
                    fout.write(line)
                    
        fout.close()
        
        # Copy the finished modified file to the environment system directory
        shutil.copy(
            fout.name,
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/step1_pdbreader.inp",
        )
        os.remove(fout.name)
        
    def _appendLigandTopparToTopparStream(self):
        # Add line to CHARMM toppar stream to read ligand parameters
        file = open(f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar.str", "a")
        file.write(f"stream {self.resname.lower()}/{self.resname.lower()}.str")
        file.close()
        
    def executeCHARMM(self, charmm_exe):
        steps = [
            "step1_pdbreader",
            "step2.1_waterbox",
            "step2.2_ions",
            "step2_solvator",
            "step3_pbcsetup_mod",
        ]

        if not self.include_ions:
            # avoid processing the step2.2_ions.inp CHARMM file
            # and remove the lines reading in these results in the
            # next CHARMM input file (step2_solvator)
            steps.remove("step2.2_ions")

            # Open a target for the modified step2_solvator input file
            fout = open(f"{self.parent_dir}/{self.ligand_id}/{self.env}/step2_solvator.inp", "wt")
            
            read_ions_section_active = False

            with open(f"{self.default_path}/step2_solvator.inp", "r+") as f:
                for line in f:
                    if line.startswith("! Add ions?"):              # Marks the start of the 
                        read_ions_section_active = True             # read ions section

                    if line.startswith("! Remove water molecules"): # Marks the end of the
                        read_ions_section_active = False            # read ions section

                    # If we are in the read ions section, we should not record these lines
                    # in the new file
                    if read_ions_section_active == True:
                        pass
                    elif line.startswith("* set"):                     # Also exclude
                        pass
                    elif line.startswith("* stream step2.2_ions.str"): # Also exclude
                        pass
                    else: # The rest of the step2_solvator.inp however, we copy verbatim
                        fout.write(line)
            fout.close()
            
            
        # Now, run CHARMM-GUI input scripts
        for step in steps:
            returncode = self._runCHARMM(step, charmm_exe)

    def _runCHARMM(self, step, charmm_exe):
        # Record the directory where we start
        main_dir = os.getcwd()
        
        # We need to go to the specific directory to run CHARMM
        os.chdir(f"{self.parent_dir}/{self.ligand_id}/{self.env}/")
        
        # Start a subprocess
        output = subprocess.run(
            [f"{charmm_exe}"] + ["-i"] + [f"{step}.inp"] + ["-o"] + [f"{step}.out"],
            text=True,
            capture_output=True,
        )
        
        # Return to the main directory
        os.chdir(main_dir)
        
        # A basic evaluation of the subprocess output codes
        if output.returncode:
            print(
                f"Something went wrong, please check the outputfile in {self.parent_dir}/{self.ligand_id}/{self.env}/{step}.out"
            )
            sys.exit("Processing terminated")
        else:
            print(f"CHARMM process finished for {step}")
            
            
    def createOpenMMSystem(self):
        print("Creating a stand-alone OpenMM system")

        # Copy files for OpenMM
        # OpenMM Python scripts
        for file in glob.glob(f"{self.default_path}/[!checkfft.py]*py"):
            shutil.copy(file, f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/")

        # Equilibration input file
        shutil.copy(
            f"{self.default_path}/omm_step4_equilibration.inp",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step4_equilibration.inp",
        )

        # Production input file
        shutil.copy(
            f"{self.default_path}/omm_step5_production.inp",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step5_production.inp",
        )

        # Toppar directory from CHARMM base directory to the OpenMM directory
        try:
            shutil.copytree(
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar",
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/toppar",
            )
        except:
            print(f"Toppar directory is already available in the OpenMM directory")

        # Convert CHARMM toppar.str to OpenMM toppar.str
        self._convertCharmmTopparStreamToOpenmm(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/toppar.str",  #
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/toppar.str",
        )

        # Inspect OpenMM toppar.str to identify ligands whose parameter
        # folders should be copied
        external_toppars = self._getExternalTopparFromTopparStream(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/toppar.str"
        )

        if external_toppars != []:
            for ext_top in external_toppars:
                try:
                    shutil.copytree(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{ext_top}",
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/{ext_top}",
                    )
                except FileNotFoundError:
                    print(
                        f"ERROR: External topology/parameter folder {ext_top} not found!"
                    )
                    
                # NOTE Here we are not allowing the overwriting of the ligand toppar folder
                # though that could be a useful behaviour. Consider adding
                # shutil.copytree(from, to, dirs_exist_ok=True) which works since Python 3.8
                except FileExistsError:
                    print(
                        f"The topology/parameters folder {ext_top} is already present."
                    )

        # Create sysinfo.dat file, which contains information about the box
        file = open(f"{self.parent_dir}/{self.ligand_id}/{self.env}/step3_pbcsetup.str")
        for line in file:
            if line.startswith(f" SET A "):
                a = line.split(" ")[-1].strip()
            elif line.startswith(f" SET B "):
                b = line.split(" ")[-1].strip()
            elif line.startswith(f" SET C "):
                c = line.split(" ")[-1].strip()
            elif line.startswith(f" SET ALPHA "):
                alpha = line.split(" ")[-1].strip()
            elif line.startswith(f" SET BETA "):
                beta = line.split(" ")[-1].strip()
            elif line.startswith(f" SET GAMMA "):
                gamma = line.split(" ")[-1].strip()

        with open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/sysinfo.dat", "w"
        ) as f:
            f.write(
                '{"dimensions":[' + f"{a}, {b}, {c}, {alpha}, {beta}, {gamma}" + "]}"
            )
        
    def _convertCharmmTopparStreamToOpenmm(self, charmmstr_path, openmmstr_path):
        new_toppar = ""
        # Read the CHARMM toppar stream from the provided location
        with open(charmmstr_path, "r") as topparstream:
            for line in topparstream:
                if line.startswith("*"):       # DROP: Header
                    pass
                elif line.startswith("!"):     # DROP: Comments
                    pass
                elif line.strip(" ") == "\n":  # DROP: Empty lines
                    pass
                elif line.startswith("read"):  # DROP: Legacy read statements
                    pass
                else:
                    strpath = line.split(" ")[-1]
                    new_toppar += f"./{strpath}"
                    # Split on whitespace and take the last item in the
                    # resulting list of words (the path to each stream file,
                    # including the newline characters.

        # Write the new OpenMM toppar stream to the requested location
        with open(openmmstr_path, "w") as topparstream:
            topparstream.write(new_toppar)

    def _getExternalTopparFromTopparStream(self, topparstream):
        external_toppars = []
        with open(topparstream, "r") as tp_str:
            for line in tp_str:
                if "toppar" in line:
                    pass
                else:
                    external_toppars.append(line.split("/")[1])

        return external_toppars


    def applyHMR(self):

        # If input_orig does not exist, the system has not been updated/overwritten yet, and
        # this is likely the first pass. Thus, we should make a backup.
        if not os.path.isfile(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input_orig.psf"
        ):

            # Copy the original PSF to a backup (input_orig)
            shutil.copy(
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input.psf",
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input_orig.psf",
            )
            print("Backed up OpenMM PSF file")

            input_psf = "step3_input.psf"

        # A backup of the orginal system exists (input_orig)
        else:
            input_psf = "step3_input_orig.psf"
            print("Reading PSF file from backup")

        ## Load parameter into ParmEd
        # PARAMETERS ARE LOADED DUE TO A BUG IN PARMED*
        #
        # This way of reading the parameters preserves the order of the
        # original/CHARMM-GUI-derived toppar.str, and ensures that
        # duplicate parameters are read in their "intended" order - the
        # last parameters read will overwrite older parameters and be
        # those used.
        #
        parms = ()
        with open(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/toppar.str", "r"
        ) as ommtopparstream:
            for line in ommtopparstream:
                parms += (line.strip("\n"),)

        # Change directory to for ParmEd to find the parameter files
        # load the parameters and change back to the original run directory
        cur_dir = os.getcwd()
        os.chdir(f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/")
        params = pm.charmm.CharmmParameterSet(*parms)
        os.chdir(cur_dir)

        # Load the PSF into ParmEd
        psf = pm.charmm.CharmmPsfFile(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/{input_psf}"
        )

        # *Since ParmEd forgot about the masses when loading the PSF above,
        # we need to load the parameters now to update the masses in the CharmmPsfFile class.
        psf.load_parameters(params)

        # Apply HMR and save PSF
        pm.tools.actions.HMassRepartition(psf).execute()
        psf.save(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/step3_input.psf",
            overwrite=True,
        )

        # Assure that hydrogen mass in a TIP3P water is greater than one
        for atom in psf:
            if atom.name.startswith("H") and atom.residue.name != "TIP3":
                assert atom.mass > 1.5
        print("Successfully applied HMR.")



    def createTFYamlFile(self, dt=0.001, nstep=2500000):

        try:
            os.makedirs(f"{self.parent_dir}/config")
        except:
            pass

        if self.env == "waterbox":

            # Creates a yaml file for ASFE simulations using TF

            fin = open(f"{self.default_path}/../temp_asfe.yaml", "r")
            fout = open(f"{self.parent_dir}/config/{self.ligand_id}.yaml", "w")
            for line in fin:
                if line.strip().startswith("name"):
                    new_line = line.replace("NAME", f"{self.ligand_id}")
                    fout.write(new_line)
                elif line.strip().startswith("tlc"):
                    new_line = line.replace("UNK", f"{self.resname}")
                    fout.write(new_line)
                elif line.strip().startswith("dt:"):
                    new_line = line.replace("TIMESTEP", str(dt))
                    fout.write(new_line)
                elif line.strip().startswith("nstep:"):
                    new_line = line.replace("NUMBEROFSTEPS", str(nstep))
                    fout.write(new_line)
                else:
                    fout.write(line)

            fin.close()
            fout.close()

            if dt != 0.001:
                self._apply_constraints(fout.name)
                
            print(f"Created YAML file {self.parent_dir}/config/{self.ligand_id}.yaml for ASFE simulations using Transformato")
                
    def _apply_constraints(self, fout):

        with open(fout, "r") as file:
            filedata = file.read()
        filedata = filedata.replace("cons: None", "cons: HBonds")
        with open(fout, "w") as file:
            file.write(filedata)
