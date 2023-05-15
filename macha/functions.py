#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in March 2023

@author: Johannes Karwounopoulos and Åsmund Kaupang
"""

import sys
import os
import glob
import shutil
import subprocess
import parmed as pm
from openbabel import openbabel
from openbabel import pybel
import string
import re
import logging
import warnings

# Supress ParmEd warnings
warnings.filterwarnings("ignore", module="parmed")


def checkInput(
    parent_dir="data", original_dir="original", protein_name=None, input_ext="pdb"
):

    # Make sure the input directory exists
    if not os.path.isdir(f"{parent_dir}/{original_dir}"):
        sys.exit(f"Input directory: {parent_dir}/{original_dir} not found!")

    # Look for an input protein for ligand exchange
    if protein_name != None:
        protein_path = f"{parent_dir}/{original_dir}/{protein_name}.{input_ext}"
        if not os.path.isfile(protein_path):
            logging.info(f"No protein for ligand exchange found.")
            logging.info(f"")
            protein_id = None
        else:
            logging.info(f"Protein for ligand exchange found at: {protein_path}")
            protein_id = os.path.splitext(os.path.basename(protein_path))[0]
            logging.info(f"Protein ID: {protein_id}")
            logging.info(f"")

    else:
        protein_id = None

    if protein_id == None:
        logging.info('No protein was specified (or found) for ligand exchange with multiple ligands.') 
        logging.info('The input is thus assumed to consist either of complexes (or double-stranded RNA)')
        logging.info('from which to create waterboxes + complexes, OR of single ligands (or single-stranded')
        logging.info('RNA) from which to create waterboxes. If your input consists of complexes, and you')
        logging.info('only wish to create waterboxes with the complexed ligands, you should use the -nc')
        logging.info('(--nocomplex) switch for efficiency.')
        logging.info(f"")
       
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
        cgenff_path,
        protein_id=None,
        rna=False,
        raw_pdb_to_mol2=False, # there is use for this as OpenBabel fails (silently) to treat certain structures e.g. triazoles, and creates a different molecule
        ligand_input_sanitation=True,
        system_ph=7.4,
        relax_input_segid_filter=False,
        #include_xray_water=False,
        segid_filter=[
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
    ):
        """
        This class prepares everything for further use with CHARMM. The PDB
        files are sliced into pieces and the ligand is converted to a mol2 file.
        A local version of CGenFF creates a stream file for the ligand.
        """
        # Instance variables
        self.parent_dir = parent_dir
        self.original_dir = original_dir
        self.ligand_id = ligand_id
        self.env: str = env
        self.cgenff_path: str = cgenff_path
        self.protein_id = protein_id
        self.rna: bool = rna
        self.raw_pdb_to_mol2: bool = raw_pdb_to_mol2
        self.ligand_input_sanitation: bool = ligand_input_sanitation
        self.system_ph = system_ph 
        self.relax_input_segid_filter: bool = relax_input_segid_filter
        #self.include_xray_water: bool = include_xray_water       
        # Other instance variables
        self.resname = str
        self.het_resnames: list = []
     
        # Class variables  
        Preparation.segidfilter = segid_filter
        
        # Relax the input segid filter if requested
        # by allowing all PRO* and HET* segids, while still barring
        # WAT*, SOLV, IONS and XRD*
        if self.relax_input_segid_filter == True:
            segfilter = Preparation.segidfilter 
            segfilter = [seg for seg in segfilter if not seg.startswith('PRO')]
            segfilter = [seg for seg in segfilter if not seg.startswith('HET')]
            Preparation.segidfilter = segfilter
            
        # WAT can also be included
        # if self.include_xray_water == True:
        #     segfilter = Preparation.segidfilter
        #     segfilter = [seg for seg in segfilter if not seg.startswith('WAT')]
        #     Preparation.segidfilter = segfilter
        # A chain-based filter will also be implemented in the future to control this function better
       
        # Make folder structure for this environment
        self.makeTFFolderStructure()
       
        # A general distinction is made base on whether a separate protein has been
        # indicated by the user
        if self.protein_id == None:
            # USE CASES
            # ligand.pdb  -> waterbox
            # complex.pdb -> complex + waterbox
            # complex.pdb -> waterbox (with -nc switch)
            
            # There is NO separate protein input file
            self.orig_protein_input = None
               
        else:
            # USE CASES
            # protein.pdb + ligand.pdb -> waterbox
            # protein.pdb + ligand.pdb -> complex + waterbox
            
            # There IS a separate protein input file
            self.orig_protein_input = (
                f"{self.parent_dir}/{self.original_dir}/{self.protein_id}.pdb"
            )
            
        # IF WE GOT HERE (CLASS INSTANTIATION) THERE IS A LIGAND/COMPLEX INPUT
        # Ligand/Complex input file
        self.orig_ligand_input = (
            f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb"
        )
        
        # Parse input
        self.inputParser()
            

        ############################################################################        
    
    def makeTFFolderStructure(self):
        self._make_folder(f"{self.parent_dir}/config")
        self._make_folder(f"{self.parent_dir}/{self.ligand_id}")
        self._make_folder(f"{self.parent_dir}/{self.ligand_id}/{self.env}")
        self._make_folder(f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm")
        self._make_folder(
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/openmm/restraints"
        )
        
    def _make_folder(self, path, silent=False):
        try:
            os.makedirs(path)
            if silent == False:
                logging.info(f"Creating folder {path} for the {self.env}")
        except OSError:
            if silent == False:
                logging.info(f"The folder {path} for the {self.env} exists - we will use it")
            if not os.path.isdir(path):
                raise
           
    def sanitizeLigandInput(self, structure_file_path):
        # This only works for PDBs
        assert structure_file_path.endswith(".pdb")
        
        # Backup original file if a backup does not exist
        if os.path.isfile(f"{structure_file_path}.org"):
            shutil.copy(
                f"{structure_file_path}.org",
                structure_file_path,
            )
        else:
            shutil.copy(
                structure_file_path,
                f"{structure_file_path}.org",
            )
        
        logging.info(f"Performing ligand input sanitation on {os.path.basename(structure_file_path)}")
        
        # Compile regular expressions
        double_cap_regex = re.compile('([A-Z][A-Z])')
                        
        # Define private functions
        def _convert_double_uppercase_to_upperlower(match_obj):
            if match_obj.group(1) is not None:
                return match_obj.group(1)[0].upper() + match_obj.group(1)[1].lower()
        
        
        # Open input and output, and write eventual changes
        # This checks and fixes:
        # - double uppercase atom names
        # - instances of duplicate atom name+number
        #
        with open(f"{structure_file_path}.org", "r") as isf:
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
                        if re.match(double_cap_regex, atomname) is not None:
                            double_cap_an = atomname.replace(' ','')
                            atomname = re.sub(double_cap_regex, _convert_double_uppercase_to_upperlower, atomname)
                            single_cap_an = atomname.replace(' ','')
                            logging.info(
                                f"Double uppercase atom name {double_cap_an} "\
                                f"replaced with {single_cap_an}"
                            )
                        
                        if re.match(double_cap_regex, elementname) is not None:
                            double_cap_en = elementname.replace(' ','').strip('\n')
                            elementname = re.sub(double_cap_regex, _convert_double_uppercase_to_upperlower, elementname)
                            single_cap_en = elementname.replace(' ','').strip('\n')
                            logging.info(
                                f"Double uppercase element name {double_cap_en} "\
                                f"replaced with {single_cap_en}"
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
                            logging.info(f"Atom name {exist_an} exists. Replaced with {new_aname + new_anumber}")

                        else:
                            # Add unique atomname to a list for comparison
                            atom_names.append(this_an)

                        # Write the finished line
                        osf.write(linestart + atomname + linemiddle + elementname)
                                                    
                    else:
                        osf.write(line)


################################################################################
################################################################################

    def inputParser(self):
        
        # USE CASES
        #
        # TYPICAL?  INPUT
        # YES       ligand.pdb  ->           waterbox (with ligand from ligand.pdb)
        # YES       complex.pdb ->           waterbox (with ligand from complex.pdb)
        # YES       complex.pdb -> complex + waterbox (with ligand from complex.pdb)
        #
        # YES       protein.pdb + ligand.pdb ->       complex + waterbox (with ligand from ligand.pdb)
        # YES       complex.pdb + ligand.pdb -> (new) complex + waterbox (with ligand from ligand.pdb)
        # NO        protein.pdb + ligand.pdb ->                 waterbox (with ligand from ligand.pdb)
        # YES       complex.pdb + ligand.pdb ->                 waterbox (with ligand from ligand.pdb)
        #

################################################################################

        # A SEPARATELY DEFINED PROTEIN IS ONLY RELEVANT IN ENVIRONMENT COMPLEX
        # SO WE WILL NOT EVEN LOOK FOR IT WHEN IN THE WATERBOX ENVIRONMENT
        # Create a list to hold the separately defined ParmEd protein objects
        sep_protein_pm_objs = []
        if self.env == "complex":
            if self.orig_protein_input != None:
                self.pm_obj = pm.load_file(
                    self.orig_protein_input, structure=True
                )
                               
                # REMOVE LONE PAIRS
                self._remove_lp()
                
                # ADD SEGMENT IDENTIFIERS
                logging.info(f"Adding segment IDs to input structure")
                self._add_segids()
                
                # Update the dataframe - IMPORTANT!
                pm_obj_df = self.pm_obj.to_dataframe()
                
                # FILTER FILTER FILTER
                # Update the object based on the input filter and update the dataframe
                # NOTE the tilde!
                self.pm_obj = self.pm_obj[~pm_obj_df.segid.isin(Preparation.segidfilter)]
                excluded_segids = ', '.join(pm_obj_df[pm_obj_df.segid.isin(Preparation.segidfilter)].segid.unique())
                if excluded_segids == "":
                    pass # excluded_segids = '[segid name missing/""]'
                else:
                    logging.info(f"Excluding segids {excluded_segids} based on segid input filter")

                pm_obj_df = self.pm_obj.to_dataframe()    
            
                # SLICE OUT THE RELEVANT PARTS - IN THIS CASE THE PROTEIN (NOT HET)
                try:
                    # The separate protein input could be a complex
                    # In this case any HET segids will be discarded
                    self.pm_obj = self.pm_obj[~pm_obj_df.segid.str.contains("HET")]
                    pm_obj_df = self.pm_obj.to_dataframe()    
                except IndexError: # An IndexError is thrown if the dataframe with which
                                   # the ParmEd object is sliced above is empty.
                    logging.error(f"Attempting to slice a ParmEd object with an empty dataframe")
                    sys.exit(1)
                    
                # To handle multiple input proteins we go through the remaining object
                # by segid
                # Start by storing the full object locally and making its dataframe
                multi_protein_obj = self.pm_obj
                multi_protein_obj_df = multi_protein_obj.to_dataframe()
                for segid in multi_protein_obj_df.segid.unique():
                    self.pm_obj = multi_protein_obj[multi_protein_obj_df.segid == segid]
                
                    # Check ionizable amino acid residue names
                    self._check_ionizable()
                
                    # We then add the protein ParmEd object to the list of protein objects
                    sep_protein_pm_objs.append(self.pm_obj)
                

################################################################################

        # NO SEPARATE PROTEIN HAS BEEN DEFINED
        if self.orig_ligand_input != None:
            self.pm_obj = pm.load_file(
                self.orig_ligand_input, structure=True
            )
            # REMOVE LONE PAIRS
            self._remove_lp()

            # ADD SEGMENT IDENTIFIERS
            if self.rna == True:
                self._add_segids_rna()
            else:
                self._add_segids()
            
            # Update the dataframe - IMPORTANT!
            pm_obj_df = self.pm_obj.to_dataframe()   
            
            # FILTER FILTER FILTER
            # Update the object based on the input filter and update the dataframe
            self.pm_obj = self.pm_obj[~pm_obj_df.segid.isin(Preparation.segidfilter)]
            excluded_segids = ', '.join(pm_obj_df[pm_obj_df.segid.isin(Preparation.segidfilter)].segid.unique())
            if excluded_segids == "":
                pass # excluded_segids = '[segid name missing/""]'
            else:
                logging.info(f"Excluding segids {excluded_segids} based on segid input filter")
            
            pm_obj_df = self.pm_obj.to_dataframe()            
                            
            # Store the "original" (filtered) ParmEd object
            # - necessary since the self.pm_obj is constantly overwritten
            pm_obj_org = self.pm_obj

#-------------------------------------------------------------------------------

            # LIGAND (SINGLE LIGAND OR FROM COMPLEX)
            # Create a list to hold all the ligand ParmEd objects
            ligand_pm_objs = []
            
            # A lone ligand or a ligand from a complex is useful in both "waterbox" and 
            # "complex" environments            
            if self.env in ["waterbox","complex"]:
                if self.rna == True:
                    # Take the RNA chain A (segid RNAA) as the guest RNA
                    try:
                        ligand_pm_obj = self.pm_obj[pm_obj_df.segid.str.contains("RNAA")]
                        ligand_pm_obj_df = ligand_pm_obj.to_dataframe()
                    except IndexError:
                        logging.error(f"No RNAA segid was found")
                        sys.exit(1)
                else:
                    try:
                        # Make a new ParmEd object and dataframe discarding any non-HET entities
                        ligand_pm_obj = self.pm_obj[pm_obj_df.segid.str.contains("HET")]
                        ligand_pm_obj_df = ligand_pm_obj.to_dataframe()
                    except IndexError:
                        logging.error(f"No HET segids were found (no ligand)")
                        sys.exit(1)
                    
                if ligand_pm_obj_df.empty:
                    logging.error(f"ERROR: No ligand or guest RNA was found in the input")
                    sys.exit(1)
                else:
                    # The reason for separating the assert statements in the below code is that
                    # in RNA many residues (with different residue numbers) will have the same chain ID. 
                    # In the case of ligands, we can sometimes meet the case of more than one ligand in
                    # the same chain (with different residue numbers), but in this case _add_segids()
                    # will have separated them, so that each HET and chain ID will match. For RNA,
                    # there will be a single HET and thus the chain ID should be "unique". To avoid 
                    # using set() and np.unique(), which are both unordered, one can use pd.unique() from
                    # pandas, or list(dict.fromkeys()) which since Python 3.5 is an ordered structure.
                    # The latter is preferred to avoid an explicit import of pandas.
                    if self.rna == True:
                        chids = list(dict.fromkeys([res.chain for res in ligand_pm_obj.residues]))
                        segids = ligand_pm_obj_df.segid.unique()
                        assert len(segids) == len(chids)
                        
                        logging.info(f"The guest RNA consists of segid/chain: "\
                            f"{', '.join(['/'.join([seg,ch]) for seg,ch in zip(segids, chids)])}"
                        )                        
                    else:
                        chids = [res.chain if res.chain != "" else "-" for res in ligand_pm_obj.residues]
                        segids = ligand_pm_obj_df.segid.unique()
                        assert len(segids) == len(chids)
                        logging.info(f"The 'ligand' consists of segid/chain: "\
                            f"{', '.join(['/'.join([seg,ch]) for seg,ch in zip(segids, chids)])}"
                        )

                if self.rna == True:
                    self.pm_obj = ligand_pm_obj
                    self._create_tlc_rna()
                    ligand_pm_objs.append(self.pm_obj)
                     
                else:
                    # Since rather often ligands in x-ray cocrystal structures carry the same
                    # residue names, we can quite brutally rename the residues, giving them an appended letter
                    # corresponding to their x in HETx
                    # BUT ONLY IF THERE IS MORE THAN ONE RESIDUE
                    if len(ligand_pm_obj.residues) > 1:
                        for idx, (resname, hetx) in enumerate([(res.name, res.segid) for res in ligand_pm_obj.residues]):
                            last_letter = hetx[-1:]
                            new_resname = resname + last_letter
                            ligand_pm_obj.residues[idx].name = new_resname
                            logging.info(f"Renamed {resname} from segid {hetx} to {new_resname}")
                        
                    # Update the dataframe
                    ligand_pm_obj_df = ligand_pm_obj.to_dataframe()
                        
                    #DEBUG
                    #logging.info(f"RESIDUE/SEGMENT NAMES: {[(res.name, res.segid) for res in ligand_pm_obj.residues]}")
                    #logging.info(f"The current HETs are: {', '.join(ligand_pm_obj_df.segid.unique())}")
                    #logging.info(f"The current resnames are: {', '.join(ligand_pm_obj_df.resname.unique())}")

                    # Beginning support for multiple HETs 
                    # which relies on that ligands in the same chain with different residue 
                    # numbers have been given different x's in HETx by _add_segids()
                    
                    # Make sure that we have understood correctly which and how many ligands there are
                    assert len([res.segid for res in ligand_pm_obj.residues]) == len([res.name for res in ligand_pm_obj.residues]) == len([res.number for res in ligand_pm_obj.residues])
                    
                    for idx, (segid, resname, resnum) in enumerate(
                        zip([res.segid for res in ligand_pm_obj.residues],
                            [res.name for res in ligand_pm_obj.residues],
                            [res.number for res in ligand_pm_obj.residues],
                        )
                    ):
                        logging.info(f"Processing {resname} {resnum} from {segid}")
                        # Add the residue name to the master list
                        self.het_resnames.append(resname)
                        
                        # Store the current segid/ligand as its own object
                        current_ligand_obj = ligand_pm_obj[ligand_pm_obj_df.segid == segid]
                        
                        # Optional ligand input (HET) sanitation
                        if ((self.ligand_input_sanitation == True) and (self.rna==False)):
                            # Save the ligand to disk, since the sanitation function works directly
                            # on the file
                            current_ligand_path = f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.ligand_id}_{segid.lower()}.pdb"
                            current_ligand_obj.save(
                                current_ligand_path, overwrite=True,
                                )

                            # Run the sanitizer, which will back up the incoming file to *.pdb.org and 
                            # make a sanitized ligand with the same name as the input)
                            self.sanitizeLigandInput(current_ligand_path)        

                            # Reload the sanitized ligand
                            sanitized_ligand_pm_obj = pm.load_file(
                                current_ligand_path, structure=True
                            )

                            # Reinstall the segid and resname that was lost in the above
                            # We are fairly sure this HET only consists of one residue, but just to be sure;
                            for res in sanitized_ligand_pm_obj.residues:
                                # Debug
                                #logging.info(f"Setting segid to {segid}")
                                res.segid = segid
                                # Debug
                                #logging.info(f"Setting resname to {resname}")
                                res.name = resname
                                # Debug
                                #logging.info(f"Setting resnum to {resnum}")                            
                                res.number = resnum
                                
                            # REINSTATE THE SANTITIZED LIGAND AS THE CENTRAL PARMED OBJECT
                            self.pm_obj = sanitized_ligand_pm_obj      
                                                      
                            # The following are common steps, regardless of ligand sanitation
                            self._remove_lp()
                            # Check for the presence of hydrogens and add if there are NONE
                            self._check_for_hydrogens() 
                            
                        else:
                            # Quick check for hydrogens if sanitation is not requested
                            if len([atm for atm in current_ligand_obj.atoms if atm.element_name == 'H']) == 0:
                                logging.warning(f"No hydrogens were found in ligand and ligand sanitation was not requested")
                                logging.warning(f"CHARMM run may fail!")
                            self.pm_obj = current_ligand_obj
                        
                        #DEBUG
                        # logging.info(f"DF AFTER LIGAND SANITATION:\n"\
                        #       f"{self.pm_obj.to_dataframe()}\n")

                        # Add the (processed) ligand objects to the list
                        ligand_pm_objs.append(self.pm_obj)
#-------------------------------------------------------------------------------
            
            # PROTEIN (FROM COMPLEX)
            # Create a list to hold the protein ParmEd objects
            protein_pm_objs = []
            # The protein part of a complex input is only useful in the environment "complex"
            if self.env == "complex":
                # We reinstall the "original" (filtered) ParmEd object and recreate its dataframe
                self.pm_obj = pm_obj_org
                pm_obj_df = self.pm_obj.to_dataframe()
                
                if self.rna == True:
                    try:
                        # Take the RNA chain B (segid RNAB) as the host RNA
                        self.pm_obj = self.pm_obj[pm_obj_df.segid.str.contains("RNAB")]
                        pm_obj_df = self.pm_obj.to_dataframe()
                    except IndexError: # An IndexError is thrown if the dataframe with which
                                    # the ParmEd object is sliced above is empty.
                        logging.error(f"Attempting to slice a ParmEd object with an empty dataframe")
                        sys.exit(1)
                        
                    if pm_obj_df.empty:
                        logging.info(f"No host RNA was found in the input from which to make a complex")
                        apo_protein_pm_obj = None
                    else:
                        assert len(pm_obj_df.segid.unique()) == len(pm_obj_df.chain.unique())

                        self._create_tlc_rna()
                        pm_obj_df = self.pm_obj.to_dataframe()
                        logging.info(f"The host RNA consists of segid/chain: "\
                            f"{', '.join(['/'.join([seg,ch]) for seg,ch in zip(pm_obj_df.segid.unique(), pm_obj_df.chain.unique())])}"
                        )
                        
                        # We then add the protein ParmEd object to the list of protein objects
                        protein_pm_objs.append(self.pm_obj)

                else:
                    try:
                        # Discard any HET segids from the ParmEd object
                        self.pm_obj = self.pm_obj[~pm_obj_df.segid.str.contains("HET")] # THROWS ERRORS IN TERMINAL IF EMPTY
                        pm_obj_df = self.pm_obj.to_dataframe()
                    except IndexError: # An IndexError is thrown if the dataframe with which
                                    # the ParmEd object is sliced above is empty.
                        logging.error(f"Attempting to slice a ParmEd object with an empty dataframe")
                        sys.exit(1)
                        
                    if pm_obj_df.empty:
                        logging.info(f"No protein was found in the input from which to make a complex")
                        apo_protein_pm_obj = None
                    else:
                        assert len(pm_obj_df.segid.unique()) == len(pm_obj_df.chain.unique())

                        logging.info(f"The 'protein' consists of segid/chain: "\
                            f"{', '.join(['/'.join([seg,ch]) for seg,ch in zip(pm_obj_df.segid.unique(), pm_obj_df.chain.unique())])}"
                        )
                        
                        # To handle multiple input proteins we go through the remaining object
                        # by segid
                        # Start by storing the full object locally and making its dataframe
                        multi_protein_obj = self.pm_obj
                        multi_protein_obj_df = multi_protein_obj.to_dataframe()
                        for segid in multi_protein_obj_df.segid.unique():
                            self.pm_obj = multi_protein_obj[multi_protein_obj_df.segid == segid]
                        
                            # Check ionizable amino acid residue names
                            self._check_ionizable()
                        
                            # We then add the protein ParmEd object to the list of protein objects
                            protein_pm_objs.append(self.pm_obj)
                        
################################################################################

        # NOW COMBINE THE PARMED OBJECT LISTS:
        self.pm_objs = []
        if self.orig_protein_input == None:
            # There is NO separately defined protein
            if self.env == "waterbox":
                self.pm_objs = self.pm_objs + ligand_pm_objs # LIST EXTENSION!
            elif self.env == "complex":
                if protein_pm_objs == []:
                    logging.error(f"")
                    logging.error(f"No protein or host RNA was found from which to make a complex")
                    logging.error(f"If you only want to make a waterbox from a complexed ligand")
                    logging.error(f"or single-stranded RNA, use the -nc switch. Mixed batches of")
                    logging.error(f"single and complexed entities are thus not supported.")
                    sys.exit("Halting the run.")
                else:
                    self.pm_objs = protein_pm_objs + ligand_pm_objs # LIST EXTENSION!
            else:
                sys.exit(f"Unrecognized environment: {self.env}")
        else:
            # There IS a separately defined protein
            if (self.env == "waterbox"):
                self.pm_objs = ligand_pm_objs # LIST EXTENSION!
            elif self.env == "complex":
                self.pm_objs = sep_protein_pm_objs + ligand_pm_objs # LIST EXTENSION!
            else:
                sys.exit(f"Unrecognized environment: {self.env}")
                
################################################################################
################################################################################


    def _remove_lp(self):
        logging.info(f"Removing lone pairs (if any)")
        # Will check if there are lone pairs and remove them
        # CGenFF will add them later on
        lps = []
        for atom in self.pm_obj:
            if atom.name.upper().startswith("LP"): # Check for lp, Lp, lP, LP
                logging.info(f"We will remove lonepair {atom.idx} {atom.name}, {atom}")
                lps.append(atom.idx)
                
        for i in range(len(lps)):
            self.pm_obj.strip(f"@{lps[i]+1-i}")
    
    def _check_ionizable(self):
        logging.info(f"Checking ionizable amino acid residue names")
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

        # Documentation trigger switches
        his_convert = False
        ion_convert = False
        
        for res in self.pm_obj.residues:  
            if res.name in ionized_res:
                ion_convert = True
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
                his_convert = True
                # TODO A lookup using PROPKA or another solution could be of interest here
                res.name = "HSD"  # ATTENTION!! Here we make all HIS to HSD
            
        # Report on what happened - could upgrade to granular reporting with residue number + name
        if ion_convert == True:
            logging.info(f"One or more ionized amino acid residue names were converted to CHARMM format")
        if his_convert == True:
            logging.warning(f"One or more HIS residues were renamed to HSD")
                
    
    def _check_for_hydrogens(self):
        # Get the segids and dataframe
        pm_obj_df = self.pm_obj.to_dataframe()
        segids = pm_obj_df.segid.unique()
        resnames = pm_obj_df.resname.unique()  
        resnums = pm_obj_df.resnum.unique()      
                
        # Give some information
        logging.info(f"Check {', '.join([' '.join([str(rna), str(rnu)]) for rna,rnu in zip(resnames,resnums)])} "\
              f"from segid {', '.join(segids)} for hydrogens"
              ) # In case there are more than one
                # resname/resnum/segid, we will see it

        # The real work starts
        if len(segids) == 1:
            segid = segids[0]
            struct = self.pm_obj[pm_obj_df.segid == segid]
            
            # Make sure there is only one residue in this segid*
            assert len(struct.residues) == 1
            
            # Work on the single residue in this segid*
            res = struct.residues[0]
            
            # Check for the presence of hydrogens in the structure                 
            if len([atm for atm in res.atoms if atm.element_name == 'H']) == 0:
                logging.info(f"No hydrogens found in residue {res.name} from segid {segid}")
                
                # Slice a ParmEd object based on the residue name (in the current segment) 
                # that does not contain any hydrogens (in case there are more than 1 residue)
                res_obj = self.pm_obj[(pm_obj_df.segid == segid) & (pm_obj_df.resname == res.name)]
                
                # Set self.resname
                self.resname = res.name
                # Debug
                #logging.info(f"THE RECORDED RESNAME IS: {self.resname}")
                
                # Make a folder for the ligand toppar files
                self._make_folder(
                    f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}",
                    silent=True
                )

                # Save the object to a PDB file
                res_obj.save(
                    f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                    f"{self.resname.lower()}/{self.resname.lower()}_noH_PE.pdb",
                    overwrite=True,
                    )
                
                # Define an OpenBabel object
                obConversion = openbabel.OBConversion()
                mol = openbabel.OBMol()
                                
                # Read a ligand PDB file to OpenBabel object mol
                obConversion.ReadFile(
                        mol,
                    f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                    f"{self.resname.lower()}/{self.resname.lower()}_noH_PE.pdb",
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
                    f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                    f"{self.resname.lower()}/{self.resname.lower()}_withH_OB.pdb",
                )
                
                # TAKE THE FIRST HALF OF THE NEW FILE FROM THE OLD (THE HEAVY ATOMS)      
                with open(
                    f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                    f"{self.resname.lower()}/{self.resname.lower()}_noH_PE.pdb",
                    "r"
                ) as ipdb_noh:                         
                    with open(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                        f"{self.resname.lower()}/{self.resname.lower()}_withH_FIX.pdb",
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
                    f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                    f"{self.resname.lower()}/{self.resname.lower()}_withH_FIX.pdb",
                    "a"
                ) as opdb:       
                    with open(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                    f"{self.resname.lower()}/{self.resname.lower()}_withH_OB.pdb",
                        "r"
                    ) as ipdb_h:
                        h_num = 1
                        resname_org = res.name
                        res.name = res.name[:3]
                        
                        for line in ipdb_h:
                            if (((line.startswith("HETATM")) or (line.startswith("ATOM"))) and (line[13:16].replace(" ","").startswith("H"))):
                                #logging.info(f"B*{line[13:25]}*{h_num}")
                                line = line.replace("UNL", res.name)
                                line = line.replace(f" H   {res.name}", f" H{h_num:<3}{res.name}")
                                line = line.replace(f"{res.name}  ", f"{res.name} {res.chain}")
                                #logging.info(f"A*{line[13:25]}*{h_num}")
                                opdb.write(line)
                                h_num += 1
                                
                        opdb.write(f"END\n")
                        res.name = resname_org
                        
                
                # Remake the ParmEd object
                self.pm_obj = pm.load_file(
                    f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                    f"{self.resname.lower()}/{self.resname.lower()}_withH_FIX.pdb",
                    structure=True
                )
                
                # Reinstall segid and resname that was lost in the above
                # We are fairly sure this HET only consists of one residue, but just to be sure;
                for res in self.pm_obj.residues:
                    # Debug
                    #logging.info(f"Setting segid to {segid}")
                    res.segid = segid
                    # Debug
                    #logging.info(f"Setting resname to {self.resname}")
                    res.name = self.resname
                    res.number = resnums[0] # Dirty (can remove after method has been secured for single ligands)
                    
                
                logging.info(
                    f"Added "\
                    f"{len([atm for atm in self.pm_obj.atoms if atm.element_name == 'H'])}"\
                    f" hydrogens to residue {res.name} for pH {self.system_ph:.1f}"
                )
                # DEBUG
                # logging.info(f"DF FROM HYDROGENATION:\n"\
                #       f"{self.pm_obj.to_dataframe()}\n")
                
                # DEBUG
                # # Resave the ParmEd object to view how it was perceived
                # self.pm_obj.save(
                #     f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                #     f"{self.resname.lower()}/{self.resname.lower()}_withH_FIX_PE.pdb",
                #     overwrite=True,
                # )
                os.remove(f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                          f"{self.resname.lower()}/{self.resname.lower()}_noH_PE.pdb"                    
                )
                os.remove(f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                          f"{self.resname.lower()}/{self.resname.lower()}_withH_OB.pdb"                    
                )
                os.remove(f"{self.parent_dir}/{self.ligand_id}/{self.env}/"\
                          f"{self.resname.lower()}/{self.resname.lower()}_withH_FIX.pdb"                    
                )

            else:
                logging.info(f"Hydrogens found in residue {res.name} from segid {segid}")
        
        else:
            sys.exit("More than 1 HET segid found!")


    def _add_segids_rna(self):
        pm_obj_df = self.pm_obj.to_dataframe()
        chids = pm_obj_df.chain.unique()
        logging.info(f"Found chain IDs: {', '.join(chids)}")
        for chain in chids:
            for res in self.pm_obj.view[
                pm_obj_df.chain == f"{chain}"
            ].residues:  # produce a select view of this chain's residues using a boolean mask
                res.segid = f"RNA{chain}"
    
    def _create_tlc_rna(self):
        # This renames the PDB residue names from single-letter codes
        # to three-letter codes ("tlc") so that ParmEd/CHARMM will understand them
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

    def _add_segids(self):
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
        
        # Get the dataframe for the current ParmEd object
        pm_obj_df = self.pm_obj.to_dataframe()

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

        rna_res = [
            "GUA",
            "CYT",
            "ADE",
            "URA",
            "THY",
            "INO",
        ]
        
        dna_res = [ # CHANGE TO DNA NAMES
            "",
        ]

        wat_res = [
            "HOH",
            "TIP",
            "SPC", # SPC waters?
        ]
        xrd_res = [# TODO Many more common crystallization additives should be added to this list
            "1PE", 
            "BOG"
        ]

        # ions_res =[# TODO A list of ions could be added for their particular handling 
        #] 

        # RUDIMENTARY TEST FOR EXISTING SEGIDS
        # The length of the ParmEd object dataframe in which the
        # segid is not blank is compared to the length of the full dataframe
      
        if len(pm_obj_df) == len(pm_obj_df[pm_obj_df.segid != ""]):
        #if len(pm_obj_df) == len(pm_obj_df.where(pm_obj_df.segid == "").dropna()): # FASTER?
            logging.info(f"All entries appear to have a segid") 
                       
        else:
            # Find the chain ids that are lacking segment ids
            chids_no_segids = pm_obj_df.where(pm_obj_df.segid == "").chain.unique()
            
            # If no chain ids are found, set chain ids based on residue names 
            # This can get messy if there are many residues - e.g. if a protein made its way here without a chain id
            # so a check is in place to limit the number upwards to 20 and stop the program.
            # PERHAPS THIS SHOULD BE MOVED TO INPUT SANITATION TO NEVER ALLOW CHAIN ID-LESS STRUCTURES IN THE DOOR
            if ((len(chids_no_segids) == 1) and (chids_no_segids == "")):
                logging.warning('No chain ID found for any of the entities that are missing segids.')
                logging.warning('Setting them to A-Z based on unique residue names (if any)')
                chid_letters = list(string.ascii_uppercase)
                resnames = pm_obj_df.resname.unique()
                if len(resnames) > 20:
                    sys.exit(f"In trying to add chain IDs to entities with missing segment IDs, based on their residue names,\n"\
                             f"more than 20 residue names were found.\n"\
                             f"Please check that a protein structure has not made its way into this method")
                    
                chids = [chid_letters.pop(0) for resn in range(len(resnames))]
                for res in self.pm_obj.view[(pm_obj_df.chain == "") & (pm_obj_df.segid == "")].residues:    
                    if res.name in resnames:
                        which_resname_idx = [idx for idx, resn in enumerate(resnames) if resn == res.name][0]
                        res.chain = chids[which_resname_idx]

                pm_obj_df = self.pm_obj.to_dataframe()
                chids_no_segids = pm_obj_df.where(pm_obj_df.segid == "").chain.unique()
                    
            # Report the findings on missing segids
            logging.info(f"Entities from chain {' and '.join(chids_no_segids)} are missing segids")
            logging.info(f"Will add segids based on residue number and chain ID")

            # Make a list of ABCD... for use as x in HETx
            # as well as a list in which to store the used "x"s
            het_letters = list(string.ascii_uppercase)
            used_het_letters = []
                
            chids = pm_obj_df.chain.unique()
           
            # FIRST PASS - FIND PROTEIN, SOLVENT AND CRYSTALLIZATION COMPONENTS
            for (
                chain
            ) in chids:
                # DEBUG
                #logging.info(f"Working on chain {chain}")

                # Loop over the residues in the current chain
                for res in self.pm_obj.view[
                    pm_obj_df.chain == f"{chain}"
                ].residues:  # produce a select view of this chain's residues using a boolean mask
                             # We use a structure view here, so that the underlying object will in
                             # fact be changed, and remember the changes
                    
                    # Separate out residues that belong to segids PRO, WAT and XRD (crystallization buffer components)
                    if res.name in aa_res:
                        res.segid = f"PRO{chain}"
                    elif res.name in wat_res:
                        res.segid = f"WAT{chain}"
                    elif res.name in xrd_res:
                        res.segid = f"XRD{chain}"
                    elif res.name in rna_res:
                        res.segid = f"RNA{chain}"
                    elif res.name in dna_res:
                        res.segid = f"DNA{chain}"
                        
            # SECOND PASS - FIND HETs
            # Update dataframe
            pm_obj_df = self.pm_obj.to_dataframe()
            
            for (
                chain
            ) in chids:

                # Debug                
                # residues_in_this_chain = pm_obj_df.where((pm_obj_df.chain == chain) & (pm_obj_df.segid == "")).dropna().resnum.unique()
                # logging.info(f"RESIDUES IN {chain}: {residues_in_this_chain}")   
                             
                # Set a temporary previous residue number
                prev_resnum = -1
                
                # Loop over a select view of this chain's residues using a boolean mask
                # We use a structure view here, so that the underlying object will in
                # fact be changed, and will remember these changes
                for res in self.pm_obj.view[(pm_obj_df.chain == chain) & (pm_obj_df.segid == "")].residues:  
                    # Note the current residue number
                    resnum = res.number
                    
                    # Remaining residues are candidates for the HET segid
                    #else:
                    
                    # There is no need to set the segid of a residue multiple times
                    # so let's check to see if we are still on the same residue
                    if resnum != prev_resnum:
                        # We primarily want to give the likely candidates
                        # for HETx an x letter that reflects, in order of priority:
                        # - a lower residue number
                        # The reason for this is that most PDBs are listed with the most important 
                        # ligand first (the one in the binding pocket), while other ligands that 
                        # appear in the crystal are given higher residue numbers
                        # This also facilitates the use of a default "PROA+HETA" filter which 
                        # makes the most common use cases easier to run
                        this_het_letter = het_letters.pop(0)
                        res.segid = f"HET{this_het_letter}"
                        used_het_letters.append(this_het_letter)
                        # Debug
                        logging.info(f"Assigned segid {res.segid} to {res.name} {res.number} from chain {res.chain}")                            
                        
                        
                        # ALTERNATIVE VARIANT, THAT PRIORITIZES CHAIN ID
                        # if chain not in used_het_letters:
                        #     res.segid = f"HET{chain}"
                        #     used_het_letters.append(str(chain))
                        #     het_letters.remove(chain)
                        #     # Debug
                        #     logging.info(f"Assigned segid {res.segid} to {res.name} {res.number} from chain {res.chain}")
                        # else:
                        #     # However, if this segid (HETx) is already used, we will 
                        #     # pop a new "x" out of the letter list. We have made sure that
                        #     # this list does not contain any letters that have been used
                        #     # previously
                        #     this_het_letter = het_letters.pop(0)
                        #     res.segid = f"HET{this_het_letter}"
                        #     used_het_letters.append(this_het_letter)
                        #     # Debug
                        #     logging.info(f"Assigned segid {res.segid} to {res.name} {res.number} from chain {res.chain}")
                    else:
                        pass      
                     
                    # Update the previous residue number                             
                    prev_resnum = resnum
            
            #logging.info(f"Added segids. Current segids: {', '.join(set(self.pm_obj.to_dataframe().segid))}")

        
    def createCRDfiles(self):
   
        used_segids = []
        for pm_obj in self.pm_objs:
            pm_obj_df = pm_obj.to_dataframe()
            segid = pm_obj_df.segid.unique()
            
            if len(segid) == 1:
                segid = str(segid[0])
            else:
                logging.error(f"ERROR: More than one segid in the ParmEd object: {segid}")
                sys.exit(1)
                
            # Store in the segids to be given to CharmmManipulation
            used_segids.append(segid)

            # WATERBOX ENVIRONMENT
            if self.env == "waterbox":
                # For ligands
                if segid.startswith("HET"):
                    # DEBUG
                    #logging.info(f"Env: {self.env}, Segid: {segid}")

                    # Note the residue name for checks
                    self.resname = (
                        pm_obj[pm_obj_df.segid == f"{segid}"].residues[0].name
                    )

                    # resname should be a 3 or 4 letter code
                    assert len(self.resname) < 5

                    # Make a folder for the ligand files in the waterbox directory
                    self._make_folder(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}",
                        silent=True
                    )
                    
                    #DEBUG
                    # logging.info(f"{self.env.upper()}: DF FOR CRD CREATION FOR {segid}:\n"\
                    #         f"{pm_obj_df[pm_obj_df.segid == segid]}\n")
                    
                    # Write the CHARMM CRD file from the ParmEd object
                    pm_obj[pm_obj_df.segid == segid].save(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                        overwrite=True,
                    )
                    
                    # # If the switch is active, copy the original PDB file to the ligand folder
                    # # This seems RISKY, since this PDB file has not been vetted by checkInputType 
                    # # (e.g for hydrogen existence)
                    # if self.small_molecule:
                    #     # We copy the pdf file since there shouldn't be any changes
                    #     shutil.copy(
                    #         f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb",
                    #         f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                    #     )
                    # else:
                    #     # If the switch is not active, 
                    #     
                    # DEBUG
                    # logging.info(f"INSPECT LIGAND DF:\n"\
                    #       f"{self.pm_obj.to_dataframe()}\n")

                    # Save a PDB file of the ligand from the ParmEd object
                    # THIS IS THE PDB FILE THAT IS CONVERTED BY OPENBABEL TO A MOL2 FILE 
                    # WHICH IS USED TO GENERATE PARAMETERS WITH CGENFF
                    pm_obj[pm_obj_df.segid == segid].save(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                        overwrite=True,
                    )
                    self.getTopparFromLocalCGenFF(self.cgenff_path)
                
                # For RNA, a segid-based naming applies (instead of resname-based as above)
                elif segid.startswith("RNA"):
                    # DEBUG
                    #logging.info(f"env: {self.env}, segid: {segid}")

                    pm_obj[pm_obj_df.segid == f"{segid}"].save(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                        overwrite=True,
                    )

            # COMPLEX ENVIRONMENT
            elif self.env == "complex":
                if segid.startswith("HET"):
                    # DEBUG
                    #logging.info(f"Env: {self.env}, Segid: {segid}")

                    # TODO: Should this check be more generally applied?
                    # Note the residue name for checks
                    self.resname = (
                        pm_obj[pm_obj_df.segid == f"{segid}"].residues[0].name
                    )
                    # resname should be a 3 or 4 letters code
                    assert len(self.resname) < 5

                    # Make a folder for the ligand in the complex directory
                    self._make_folder(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}",
                        silent=True
                    )
                    
                    #DEBUG             
                    # logging.info(f"{self.env.upper()}: DF FOR CRD CREATION FOR {segid}:\n"\
                    #         f"{pm_obj_df[pm_obj_df.segid == segid]}\n")
                    
                    # Save the CHARMM CRD file of the ligand from the ParmEd object
                    pm_obj[pm_obj_df.segid == segid].save(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                        overwrite=True,
                    )
                    
                    # Save a PDB file of the ligand from the ParmEd object
                    # THIS IS THE PDB FILE THAT IS CONVERTED BY OPENBABEL TO A MOL2 FILE 
                    # WHICH IS USED TO GENERATE PARAMETERS WITH CGENFF
                    pm_obj[pm_obj_df.segid == f"{segid}"].save(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
                        overwrite=True,
                    )
                    self.getTopparFromLocalCGenFF(self.cgenff_path)

                else: 
                    # DEBUG
                    #logging.info(f"Env: {self.env}, Segid: {segid}")

                    # Save the protein, as well as eventual other segids,
                    # as CHARMM CRD and -PDB files, directly in the complex directory
                    pm_obj[pm_obj_df.segid == f"{segid}"].save(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.crd",
                        overwrite=True,
                    )
                    pm_obj[pm_obj_df.segid == f"{segid}"].save(
                        f"{self.parent_dir}/{self.ligand_id}/{self.env}/{segid.lower()}.pdb",
                        overwrite=True,
                    )
            else:
                # Alert if unsupported environments are requested or if input has typos
                sys.exit(f"Unrecognized environment: {self.env}")

        return used_segids
    
    def getTopparFromLocalCGenFF(
        self,
        cgenff_path=False,
    ):
                
        cgenff_bin = None
        cgenff_output = None

        # If no particular path is given, check whether CGenFF is available
        if ((cgenff_path == False) or (cgenff_path == "")):
            cgenff_path = shutil.which("cgenff")
            if cgenff_path == None:
                logging.error("")
                logging.error("This package requires cgenff for parameterization.")
                logging.error("Please install it in the active environment or point the routine")
                logging.error("to the right path using the key cgenff_path='/path/to/cgenff' .")
                sys.exit(1)
            else:
                cgenff_bin = cgenff_path
        else:
            cgenff_bin = cgenff_path

        # CGenFF exists - start program
        if cgenff_bin != None:
            
            # CGenFF needs a mol2 file as input file
            self._create_mol2_sdf_file()
            #self._create_mol2_sdf_file_pybel() #alternative
            
            # If a ligand toppar stream file exists (from a previous run), 
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
                logging.error(f"CGenFF returned an error after being called with:")
                logging.error(f"{' '.join(cgenff_output.args)}")
                logging.error(f"{cgenff_output.stdout}")
                logging.error(f"{cgenff_output.stderr}")
                #sys.exit(1)
            else:
                logging.info(f"CGenFF executed successfully")

            # Modify the resname in the new stream file
            self._modify_resname_in_stream()

        return cgenff_output

    def _create_mol2_sdf_file_pybel(self):
        mol = next(
            pybel.readfile(
                "pdb", 
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb"
            )
        )
        # These may be helpful, but seem to erease the atom numbers (in the MOL2 file)
        # which causes BILDC to fail having been called with a null IC table
        #mol.OBMol.ConnectTheDots()
        #mol.OBMol.PerceiveBondOrders()
        
        ob_reschains = [str(res.GetChain()) if res.GetChain() != " " else "-" for res in openbabel.OBResidueIter(mol)]
        ob_resnames  = [str(res.GetName()) if str(res.GetName()) != " " else "-" for res in openbabel.OBResidueIter(mol)]
        ob_resnumatoms = [str(res.GetNumAtoms()) if str(res.GetNumAtoms()) != " " else "-" for res in openbabel.OBResidueIter(mol)]

        logging.info(f"Converting ligand from {self.ligand_id} ({mol.GetFormula()}) to MOL2 and SDF files")
        logging.info(f"The OBMol object consists of the following residues (chain/name/#atoms): "\
                     f"{', '.join([f'{chid}/{name}/{numa}' for chid,name,numa in zip(ob_reschains, ob_resnames, ob_resnumatoms)])}")
                    
        assert (mol.OBMol.NumResidues()) == 1
        
        mol.write(
            "mol2",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.mol2"
            )
        mol.write(
            "sdf",
            f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.sdf"
            )

    def _create_mol2_sdf_file(self):

        obConversion = openbabel.OBConversion()
        mol = openbabel.OBMol()
                
        # Read a ligand PDB file to OpenBabel object mol
        
        # THIS APPEARS TO BE REDUNDANT SINCE small_molecule == True already copies 
        # the original PDB to the ligand folder
        
        # From the "original" directory
        if self.raw_pdb_to_mol2 == True:
            logging.info(f"ATTENTION: Using raw PDB file directly from input for conversion to MOL2 and toppar generation")
            logging.info(f"This requires a fully hydrogenated ligand and will fail if the input PDB is a complex")
            obConversion.ReadFile(
                mol,
                f"{self.parent_dir}/{self.original_dir}/{self.ligand_id}.pdb",
            )
        # From the "ligandID/environment/residueName" directory, e.g. "cpd10/waterbox/dy6/dy6.pdb"
        # That is, from the PDB from the ParmEd object
        else:
            obConversion.ReadFile(
                mol,
                f"{self.parent_dir}/{self.ligand_id}/{self.env}/{self.resname.lower()}/{self.resname.lower()}.pdb",
            )

        ob_reschains = [str(res.GetChain()) if res.GetChain() != " " else "-" for res in openbabel.OBResidueIter(mol)]
        ob_resnames  = [str(res.GetName()) if str(res.GetName()) != " " else "-" for res in openbabel.OBResidueIter(mol)]
        ob_resnumatoms = [str(res.GetNumAtoms()) if str(res.GetNumAtoms()) != " " else "-" for res in openbabel.OBResidueIter(mol)]

        logging.info(f"Converting ligand from {self.ligand_id} ({mol.GetFormula()}) to MOL2 and SDF files")
        logging.info(f"The OBMol object consists of the following residues (chain/name/#atoms): "\
                     f"{', '.join([f'{chid}/{name}/{numa}' for chid,name,numa in zip(ob_reschains, ob_resnames, ob_resnumatoms)])}")
                            
        assert (mol.NumResidues()) == 1
        
        # These may be helpful, but seem to erease the atom numbers (in the MOL2 file)
        # which causes BILDC to fail having been called with a null IC table
        #mol.ConnectTheDots()
        #mol.PerceiveBondOrders()
        
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

 
