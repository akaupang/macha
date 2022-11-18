#!/bin/bash

# System variables
# Macha python backend base
macha_py_base="/home/manny/Documents/Work/UiO/Modeling/wien/scripts/bash/macha/macha"

# CHARMM binary
# NOTE - if a path to the CHARMM binary
# Check if "charmm" is in the PATH
charmm_bin_man="/home/manny/charmm/bin/charmm"

if [ "${charmm_bin_man}" = "" ];then
  if ! command -v charmm &> /dev/null; then
    echo "The CHARMM binary must be available to this tool, either through"
    echo "the PATH, or by setting its location in the script."
    echo ""
    echo "Current PATH:"
    echo $PATH
    exit
  else
    charmm_bin="charmm"
  fi
else
  charmm_bin=${charmm_bin_man}
fi

# MBL system modification script directory
mblpyloc_def="/home/manny/Documents/Work/UiO/Modeling/wien/proteins/mbls/patches/python"
mblpyname_mm="mbl_system_modifications.py"
mblpyname_qmmm="mbl_qmmm_setup_v2.py"
mblconfname="mbl_mod_configuration.yaml"

# Location of CHARMM-GUI-provided OpenMM Python scripts
# Default is the directory created upon creation of a modified OpenMM system,
# that is the original "openmm" directory
cgommpyloc_def="$PWD/prev_openmm_1"

# The name of the internally generated step 3.1 CHARMM input script file is
ommcname="step3.1_omm"

# The name of the internally generated HMR application Python script (which uses ParmEd) is
ahmrf_py="apply_HMR_with_ParmEd.py"

# TRANSFORMATO_RELATED
# ALL FUNCTIONS IN PYTHON BACKEND

# Functions
# Execute CHARMM
charmm_exec () {
  ${charmm_bin} "$@"
}

# Execute CGenFF (local binary)
cgenff_exec () {
  ${cgenff_bin} "$@"
}

# Transformato submenu
run_t_submenu () {
  echo ""
  echo "Transformato Options"
  echo ""
  echo "Key"
  echo " 1 : Run Transformato system creation (Python 3.9)"
  echo " 2 : Option2"
  echo " 3 : Option3"
  echo " 4 : Option4"
  echo ""
  echo " r : Return to main menu"
  echo " q : Quit"
  echo ""
  read -n 1 -p "Process options (1 - 4, r, q): " "submenuoption"
################################################################################

  if [ "$submenuoption" = "1" ]; then
    echo ""
    echo "Please provide the base path for transformato"
    call_path=$PWD
    echo "Press Enter to use the current directory:"
    echo ${call_path}
    echo ""
    echo "or press d     to provide your own directory."
    read -r -s -n 1 key

    if [ "$key" = "" ]; then
      tf_base=$call_path
    elif [ "$key" = "d" ]; then
      echo "Enter new directory"
      read tf_base      
      if [ "$tf_base" = "" ]; then
        run_t_submenu
      else
        tf_base=$(realpath "$tf_base")
        cd $tf_base
      fi
    else
      run_t_submenu
    fi  

    # Prepare for python backend run  
    # Copy main.py from default location
    cp $macha_py_base/main.py ./

    # Edit import statement
    sed -i "s/from functions import/from macha.functions import/g" main.py

    # Run backend
    python3 main.py
    
    # Return to Transformato submenu
    run_t_submenu

  elif [ "$submenuoption" = "2" ]; then
    run_t_submenu

  elif [ "$submenuoption" = "r" ];then
    runmenu

  elif [ "$submenuoption" = "q" ];then
    echo " <- Goodbye"
    exit 0       

  else
    timeout 2s echo " <- Unrecognized option. Press any key."
    #read -n 1 k <&1
    run_t_submenu
  fi
}

# Advanced submenu
run_a_submenu () {
  echo ""
  echo "Advanced Options"
  echo ""
  echo "Key  Pre-run options"
  echo " 1 : Add ligand toppar to toppar.str"
  echo " 2 : Modify MM input files with MBL customization hooks (Python 3.9)"
  echo " 3 : Create ligand toppar folders and parameters with local CGenFF"
  echo " 4 : Create sQM/MM input file for MBLs using SCCDFTB"
  echo " 5 : Create QM/MM input file for MBLs using Q-CHEM"
  echo " 6 : Create (s)QM/MM input file for MBLs with custom input/output"
  echo ""
  echo "Key  Post-run options"
  echo " 7 : Apply HMR to OpenMM system"
  echo ""
  echo " r : Return to main menu"
  echo " q : Quit"
  echo ""
  read -n 1 -p "Process options (1 - 4, r, q): " "submenuoption"
################################################################################

  if [ "$submenuoption" = "1" ]; then
    echo ""
    echo "Input ligand name=ligand toppar directory name (or type \"return\"+Enter to return)"
    read ligname
    if [ "${ligname}" = "return" ]; then
      run_a_submenu
    else
      ligtopparblock="\n"
      ligtopparblock+="! Ligand topology and parameter files\n"
      ligtopparblock+="open read card unit 10 name ${ligname}/${ligname}.rtf\n"
      ligtopparblock+="read  rtf card unit 10 append\n"
      ligtopparblock+="\n"
      ligtopparblock+="open read card unit 20 name ${ligname}/${ligname}.prm\n"
      ligtopparblock+="read para flex card unit 20 append\n"
      ligtopparblock+="\n"

      printf "%b" "$ligtopparblock" >> "toppar.str"
      
      echo ""
      echo ""
      echo "Added ligand topology/parameter files to the toppar.str"
    fi
    run_a_submenu

  elif [ "$submenuoption" = "2" ];then
    mblpyloc_def=$(realpath "$mblpyloc_def")
    echo ""
    echo ""
    echo "Looking for MBL system modification Python scripts in:"
    echo ${mblpyloc_def}
    if [[ -f "${mblpyloc_def}/${mblpyname_mm}" && -f "$PWD/${mblconfname}" ]];then
      python3 -u ${mblpyloc_def}/${mblpyname_mm} $PWD ${mblconfname} &&\
      echo "System modification hooks added to input files"
    else
      echo "MBL system modification scripts not found! Please provide a valid location and set this in the macha script file."
    fi

    run_a_submenu

  elif [ "$submenuoption" = "3" ];then
    echo ""
    echo "Provide a path to the ligands in MOL2 format"
    read lig_path
    lig_path=$(realpath "$lig_path")
    call_path=$PWD

    cd $lig_path
    IFILES="$lig_path/*.mol2"
    for ifile in ${IFILES}; do
      if [ -f "${ifile}" ]; then
        filename=$(basename -- "${ifile}")
        name="${filename%.*}"
        rm -rf ${name} &&\   # Remove a previous directory with this name
        mkdir -p ${name} &\ # Make a new directory
        cp ${ifile} ${name}/ # Copy the ligand there
        (cd $lig_path/${name} &&\
         cgenff_exec ${filename} -v -f "${name}.str" -m "${name}.log" &&\
         echo "Processed ${name} and created a correspondingly named toppar directory")
      fi
    done
    cd call_path
    run_a_submenu

  elif [ "$submenuoption" = "4" ];then
    mblpyloc_def=$(realpath "$mblpyloc_def")
    echo ""
    echo ""
    echo "Looking for MBL system modification Python scripts in:"
    echo ${mblpyloc_def}
    if [[ -f "${mblpyloc_def}/${mblpyname_qmmm}" && -f "$PWD/${mblconfname}" ]];then
      python3 -u ${mblpyloc_def}/${mblpyname_qmmm} $PWD ${mblconfname} "sccdftb"
    else
      echo "MBL system modification Python scripts or YAML configuration file not found! Please provide valid locations and set these at the top of the macha script file."
    fi

    run_a_submenu
    
  elif [ "$submenuoption" = "5" ];then
    mblpyloc_def=$(realpath "$mblpyloc_def")
    echo ""
    echo ""
    echo "Looking for MBL system modification Python scripts in:"
    echo ${mblpyloc_def}
    if [[ -f "${mblpyloc_def}/${mblpyname_qmmm}" && -f "$PWD/${mblconfname}" ]];then
      python3 -u ${mblpyloc_def}/${mblpyname_qmmm} $PWD ${mblconfname} "qchem"
    else
      echo "MBL system modification Python scripts or YAML configuration file not found! Please provide valid locations and set these at the top of the macha script file."
    fi

    run_a_submenu
    
  elif [ "$submenuoption" = "6" ];then
    call_path=$PWD
    echo ""
    echo ""  
    echo "Provide a path to the source system folder"
    echo "Press Enter to use the current directory:"
    echo ${call_path}
    echo ""
    echo "or press d     to provide your own directory."
    read -r -s -n 1 key
    if [ "$key" = "" ]; then
      sys_base=$call_path
    elif [ "$key" = "d" ]; then
      echo "Enter new directory"
      read sys_base      
      if [ "$sys_base" = "" ]; then
        run_a_submenu
      else
        sys_base=$(realpath "$sys_base")
        cd $sys_base
      fi
    else
      run_a_submenu
    fi

    echo ""
    echo "Provide the base name of the source PSF"
    echo "Press Enter to use 'step1_pdbreader':"
    echo ""
    echo "or press d     to provide your own name."
    read -r -s -n 1 key
    if [ "$key" = "" ]; then
      psf_base="step1_pdbreader"
    elif [ "$key" = "d" ]; then
      echo "The following PSF files are available in the source directory:"
      for entry in "$sys_base"/*.psf; do
        name_ext=$(basename -- "${entry}")
        echo "${name_ext%.*}"
      done
      echo "Enter PSF base name"
      read psf_base      
      if [ "$psf_base" = "" ]; then
        run_a_submenu
      else
        psf_base="${psf_base}"
      fi
    else
      run_a_submenu
    fi
    
    echo ""
    echo "The following CRD files are available in the source directory:"
    for entry in "$sys_base"/*.crd; do
      name_ext=$(basename -- "${entry}")
      echo "${name_ext%.*}"
    done
    echo "Enter CRD base name"

    read crd_base
    if [ "$crd_base" = "" ]; then
      run_a_submenu
    else
      crd_base="${crd_base}"
    fi
    
    echo ""
    echo "Create input files for sccdftb or qchem?"
    echo "Press Enter to use 'qchem':"
    echo ""
    echo "or press d     to set it manually."
    read -r -s -n 1 key
    if [ "$key" = "" ]; then
      theorylevel="qchem"
    elif [ "$key" = "d" ]; then
      echo "Choose theory level (sccdftb or qchem)"
      read theorylevel      
      if [ "$theorylevel" = "" ]; then
        run_a_submenu
      else
        theorylevel="${theorylevel}"
      fi
    else
      run_a_submenu
    fi

    echo ""
    echo "Provide the base name for the new input files"
    echo "Press Enter to use 'step1c':"
    echo ""
    echo "or press d     to provide your own name."
    read -r -s -n 1 key
    if [ "$key" = "" ]; then
      new_base="step1c"
    elif [ "$key" = "d" ]; then
      echo "Enter new base name"
      read new_base      
      if [ "$new_base" = "" ]; then
        run_a_submenu
      else
        new_base="${new_base}"
      fi
    else
      run_a_submenu
    fi
    
    mblpyloc_def=$(realpath "$mblpyloc_def")
    echo ""    
    echo "Looking for MBL system modification Python scripts in:"
    echo ${mblpyloc_def}
    if [[ -f "${mblpyloc_def}/${mblpyname_qmmm}" && -f "${sys_base}/${mblconfname}" ]];then
      python3 -u ${mblpyloc_def}/${mblpyname_qmmm} ${sys_base} ${mblconfname} ${theorylevel} ${psf_base} ${crd_base} ${new_base}
    else
      echo "MBL system modification Python scripts or YAML configuration file not found! Please provide valid locations and set these at the top of the macha script file."
    fi

    run_a_submenu

  elif [ "$submenuoption" = "7" ];then
    
    ahmrf_str='import os                                                                   \n'
    ahmrf_str+='import sys                                                                  \n'
    ahmrf_str+='import shutil                                                               \n'
    ahmrf_str+='import parmed as pm                                                         \n'
    ahmrf_str+='import warnings                                                             \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='# Supress ParmEd warnings                                                   \n'
    ahmrf_str+='warnings.filterwarnings("ignore", module="parmed")                          \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='try:                                                                        \n'
    ahmrf_str+='    input_dir = sys.argv[1]                                                 \n'
    ahmrf_str+='except IndexError:                                                          \n'
    ahmrf_str+='    input_dir = "."                                                         \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='# If the system has not been updated/overwritten (input_orig does not exist)\n'
    ahmrf_str+='if not os.path.isfile(f"{input_dir}/openmm/step3_input_orig.psf"):          \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='    # Copy the original PSF to a backup (input_orig)                        \n'
    ahmrf_str+='    shutil.copy(f"{input_dir}/openmm/step3_input.psf",#                     \n'
    ahmrf_str+='                f"{input_dir}/openmm/step3_input_orig.psf")                 \n'
    ahmrf_str+='    print("Backed up OpenMM PSF file")                                      \n'
    ahmrf_str+='    input_psf = "step3_input.psf"                                           \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='else:                                                                       \n'
    ahmrf_str+='    input_psf = "step3_input_orig.psf"                                      \n'
    ahmrf_str+='    print("Reading PSF file from backup")                                   \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='# Read parameters from toppar.str to preserve loading order                 \n'
    ahmrf_str+='parms = ()                                                                  \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='with open(f"{input_dir}/toppar.str", "r") as topparstream:                  \n'
    ahmrf_str+='    for line in topparstream:                                               \n'
    ahmrf_str+='        if line.startswith("*"):          # Header                          \n'
    ahmrf_str+='            pass                                                            \n'
    ahmrf_str+='        elif line.startswith("!"):        # Comments                        \n'
    ahmrf_str+='            pass                                                            \n'
    ahmrf_str+='        elif line.strip(" ") == "\\n":     # Empty lines                     \n' #NOTE ESCAPE OF BACKSLASH
    ahmrf_str+='            pass                                                            \n'
    ahmrf_str+='        elif line.startswith("read"):     # Legacy read statements          \n'
    ahmrf_str+='            pass                                                            \n'
    ahmrf_str+='        else:                                                               \n'
    ahmrf_str+='            line = line.strip("\\n")                                         \n' #NOTE ESCAPE OF BACKSLASH
    ahmrf_str+='            parms += ( line.split(" ")[-1], )                               \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='params = pm.charmm.CharmmParameterSet(*parms)                               \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='# Load the PSF into ParmEd                                                  \n'
    ahmrf_str+='psf = pm.charmm.CharmmPsfFile(f"{input_dir}/openmm/{input_psf}")            \n'
    ahmrf_str+='psf.load_parameters(params)                                                 \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='# Apply HMR and save PSF                                                    \n'
    ahmrf_str+='pm.tools.actions.HMassRepartition(psf).execute()                            \n'
    ahmrf_str+='psf.save(f"{input_dir}/openmm/step3_input.psf", overwrite = True)           \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='# Assure that hydrogen masses are greater than one                          \n'
    ahmrf_str+='for atom in psf:                                                            \n'
    ahmrf_str+='    if atom.name.startswith("H") and atom.residue.name != "TIP3":           \n'
    ahmrf_str+='        assert atom.mass > 1.5                                              \n'
    ahmrf_str+='                                                                            \n'
    ahmrf_str+='print("Successfully applied HMR.")                                          \n'
      
    # Determine input directory
    echo ""
    echo "Please provide the base path for the system. The directory should contain"
    echo "toppar and openmm directories, as well as the toppar stream file."
    call_path=$PWD
    echo "Press Enter to use the current directory:"
    echo ${call_path}
    echo ""
    echo "or press d     to provide your own directory."
    read -r -s -n 1 key

    if [ "$key" = "" ]; then
      sys_base=$call_path
    elif [ "$key" = "d" ]; then
      echo "Enter new directory"
      read sys_base      
      if [ "$sys_base" = "" ]; then
        run_a_submenu
      else
        sys_base=$(realpath "$sys_base")
        cd $sys_base
      fi
    else
      run_a_submenu
    fi  

    # Write python file
    printf "%b" "${ahmrf_str}" > "${ahmrf_py}"

    # Apply HMR
    python3 ${ahmrf_py} ${sys_base}

    # Return to submenu  
    run_a_submenu

  elif [ "$submenuoption" = "r" ];then
    runmenu

  elif [ "$submenuoption" = "q" ];then
    echo " <- Goodbye"
    exit 0       

  else
    timeout 2s echo " <- Unrecognized option. Press any key."
    #read -n 1 k <&1
    run_a_submenu
  fi
}

# Main menu
runmenu () {
  echo ""
  echo "Manual Launch Control for CHARMM-GUI Scripts"
  echo ""
  echo "Key"
  echo " a : Advanced pre-/post-run options"
  echo " t : Transformato options"
  echo ""
  echo " 1 : Step 1   PDB Reader (modified)"
  echo " 2 : Step 2.1 Waterbox"
  echo " 3 : Step 2.2 Ions"
  echo " 4 : Step 2   Solvator"
  echo " 5 : Step 3   PBC Setup"
  echo ""
  echo " c : Run Steps 1-3 consecutively"
  echo ""
  echo " 6 : Step 3.1 Convert CHARMM system to a stand-alone OpenMM system"
  echo " 7 : Step 3.2 Copy CG OpenMM CHARMM interpreter python scripts to openmm/"
  echo " 8 : Step 3.3 Set system base name"
  echo "              Set simulation length (and timestep)"
  echo "              Create replicas"
  echo ""
  echo " q : Quit"
  echo ""
  read -n 1 -p "Process options (1 - 9, q): " "menuoption"

###############################################################################

  if [ "$menuoption" = "1" ]; then
    echo ""
    charmm_exec < step1_pdbreader.inp > step1_pdbreader.out && tail -n 6 step1_pdbreader.out
    sed -n '/START_PAR/,/END_PAR/p' step1_pdbreader.out > step1_used_parameters.dat
    sed -n '/START_HBUILD/,/END_HBUILD/p' step1_pdbreader.out > step1_hbuild_log.dat
    runmenu
  elif [ "$menuoption" = "2" ]; then
    echo ""
    charmm_exec < step2.1_waterbox.inp > step2.1_waterbox.out && tail -n 6 step2.1_waterbox.out && runmenu

  elif [ "$menuoption" = "3" ]; then
    echo ""
    charmm_exec < step2.2_ions.inp > step2.2_ions.out && tail -n 6 step2.2_ions.out && runmenu

  elif [ "$menuoption" = "4" ]; then
    echo ""
    charmm_exec < step2_solvator.inp > step2_solvator.out && tail -n 6 step2_solvator.out && runmenu

  elif [ "$menuoption" = "5" ]; then
    echo ""
    charmm_exec < step3_pbcsetup.inp > step3_pbcsetup.out && tail -n 6 step3_pbcsetup.out && runmenu

  elif [ "$menuoption" = "c" ]; then
    echo ""
    echo "Step 1 PDB Reader:"
    charmm_exec < step1_pdbreader.inp > step1_pdbreader.out &&\
      tail -n 6 step1_pdbreader.out &&\
    echo "Step 2.1 Waterbox:" &&\
    charmm_exec < step2.1_waterbox.inp > step2.1_waterbox.out &&\
      tail -n 6 step2.1_waterbox.out &&\
    echo "Step 2.2 Ions:" &&\
    charmm_exec < step2.2_ions.inp > step2.2_ions.out &&\
      tail -n 6 step2.2_ions.out &&\
    echo "Step 2 Solvator" &&\
    charmm_exec < step2_solvator.inp > step2_solvator.out &&\
      tail -n 6 step2_solvator.out &&\
    echo "Step 3 PBC Setup" &&\
    charmm_exec < step3_pbcsetup.inp > step3_pbcsetup.out &&\
      tail -n 6 step3_pbcsetup.out &&\
    sed -n '/START_PAR/,/END_PAR/p' step1_pdbreader.out > step1_used_parameters.dat
    runmenu

  elif [ "$menuoption" = "6" ]; then
    echo ""
    echo ""
    # Handle an (eventual) existing openmm folder, to allow backup and continued modification (rerun of step3.1). Also handle the case of a new folder (files were recopied, system reset). NOTE may be moot
    if [ -d "openmm" ]; then
      echo "Found existing OpenMM directory (openmm/)"
      
      # Look for previously backed up folders and give new backup a higher number
      podircount="$(find ./prev_openmm* -maxdepth 0 -mindepth 0 -type d -printf x 2>/dev/null | wc -c)"
      podirnum=`echo "1 + $podircount" | bc`
      podir="prev_openmm_$podirnum"
      mkdir -p "$podir/"

      mv openmm/* "$podir/"
      echo "Backed it up as $podir/"

      # Check if a template directory is present # Consider if this option is moot!
      if [ -d "new_openmm" ]; then
        `/bin/rm -rf openmm/`
        mv new_openmm/ openmm/
        mkdir -p openmm/restraints/
        echo "Created an OpenMM directory (openmm/) from template files"
      else
        mkdir -p openmm/
        mkdir -p openmm/restraints/
        # Transformato styled openMM folder
        mkdir -p openmm/toppar/
        cp -r toppar/* openmm/toppar/
        echo "Created a new OpenMM directory for the current system"
      fi
    else
      if [ -d "new_openmm" ]; then
        mv new_openmm/ openmm/
        echo "Created an OpenMM directory (openmm/) from template files"
      else
        mkdir -p openmm
        echo "Created a new OpenMM directory for the current system"
      fi
      mkdir -p openmm/restraints/
    fi

    # TOPPAR STREAM AND ADDED TOPOLOGIES/PARAMETERS
    # Transform the CHARMM-GUI toppar.str to an OpenMM toppar.str
    # at least with the formers current formatting
    # The following sed command searches for anything + "toppar/" 
    # and prints out "../" + the rest of the line on which it found 
    # the search parameter, which is now; "line ends with rtf/prm/str".
    `grep -v "!" toppar.str > toppar_tmp.str`

    #`sed 's/.* //' toppar_tmp.str | sed -n -e 's/^.*\(rtf$\|prm$\|str$\)/..\/\0/p' > openmm/toppar.str` &&\
    # Transformato styled toppar.str
    `sed 's/.* //' toppar_tmp.str | sed -n -e 's/^.*\(rtf$\|prm$\|str$\)/\0/p' > openmm/toppar.str` &&\
    `/bin/rm -f toppar_tmp.str` &&\
    echo "Converted toppar.str (CHARMM format) to openmm/toppar.str (OpenMM/Transformato format)"
    # Find lines in openmm/toppar.str that do not contain the word "toppar"
    # These are likely components with topology/parameters added separately
    # which will need to be copied to the openmm folder
    # The following command processes those lines by keeping everything up 
    # to the (last) forward slash and then filtering these by uniqueness
    # to remove duplicate entries originating from separate reading of
    # top/rtf- and par/prm files.
    add_ids=`grep -v "toppar" openmm/toppar.str | sed -e 's:^\(.*\)\(/.*\)$:\1:' | uniq`
    
    if [ ! "$add_ids" = "" ]; then
      # Set the separator temporarily to '\n' and create and array from the ids
      IFS=$'\n' add_ids_arr=(${add_ids})
    
      # Loop through the ids and copy their folders to the openmm folder
      for thisid in "${add_ids_arr[@]}"; do
        `cp -r $thisid openmm/`
      done
    else
      echo "No additional topologies/parameters found in toppar.str"
    fi

    # Copy toppar folder to the OpenMM directory if it does not already exist
    if [ -d "openmm/toppar" ]; then
      echo "Toppar folder already present."
    else
      `cp -r toppar/ openmm/toppar` &&\
      echo "Copied toppar folder to the OpenMM directory."
    fi

    # SIMULATION BOX SIZE
    # Produce sysinfo.dat from step2.1_waterbox.prm 
    # Only for cubic box so far
    dima=`sed -n -e "s/^.*\(SET A = \)//p" step2.1_waterbox.prm`
    dimb=`sed -n -e "s/^.*\(SET B = \)//p" step2.1_waterbox.prm`
    dimc=`sed -n -e "s/^.*\(SET C = \)//p" step2.1_waterbox.prm`
    anga=`sed -n -e "s/^.*\(SET ALPHA = \)//p" step2.1_waterbox.prm`
    angb=`sed -n -e "s/^.*\(SET BETA  = \)//p" step2.1_waterbox.prm`
    angc=`sed -n -e "s/^.*\(SET GAMMA = \)//p" step2.1_waterbox.prm`

    echo "{\"dimensions\": [$dima.0, $dimb.0, $dimc.0, $anga, $angb, $angc]}" > openmm/sysinfo.dat &&\
    echo "Produced sysinfo.dat from step2.1_waterbox.prm"

    # OMM Conversion CHARMM Input File
    ocif=""
    ocif+="* CONVERSION TO AN OPENMM SYSTEM                                                \n"
    ocif+="* REVERSE-ENGINEERED FROM CHARMM-GUI SCRIPTS                                    \n"
    ocif+="* ONLY CUBIC BOXES SUPPORTED!                                                   \n"
    ocif+="* SCRIPT MAY NEED MODIFICATION TO SUIT YOUR SYSTEM!                             \n"
    ocif+="                                                                                \n"
    ocif+="DIMENS CHSIZE 5000000 MAXRES 3000000                                            \n"
    ocif+="                                                                                \n"
    ocif+="! Read topology and parameter files                                             \n"
    ocif+="stream toppar.str                                                               \n"
    ocif+="                                                                                \n"
    ocif+="! Read PSF and coordinates                                                      \n"
    ocif+="read psf  card name step3_pbcsetup.psf                                          \n"
    ocif+="read coor card name step3_pbcsetup.crd                                          \n"
    ocif+="                                                                                \n"
    ocif+="! Transform coordinates to OpenMM style                                         \n"
    ocif+="coor stat                                                                       \n"
    ocif+="calc xr = -1*?XMIN                                                              \n"
    ocif+="calc yr = -1*?YMIN                                                              \n"
    ocif+="calc zr = -1*?ZMIN                                                              \n"
    ocif+="coor tran xdir @XR ydir @YR zdir @ZR                                            \n"
    ocif+="                                                                                \n"
    ocif+="! Positional restraints for equilibration                                       \n"
    ocif+="! from CG step4_equilibration.inp                                               \n"
    ocif+="! MAY NEED MODIFICATION TO COVER THE ACTUAL SYSTEM                              \n"
    ocif+="! Define BB and SC (and likely more in other systems)                           \n"
    ocif+="! BB contains:                                                                  \n"
    ocif+="! - selected heavy atoms, apparently the backbone of the protein (PROT)         \n"
    ocif+="!   notably not; OT1,OT2,                                                       \n"
    ocif+="! - selected heavy atoms in sugars (CARB) [complicated by "bonded" - pyranose?] \n"
    ocif+="! - all non-hydrogens in ligands or metals (HETE)                               \n"
    ocif+="!                                                                               \n"
    ocif+="! SC contains those heavy atoms not selected in BB, from PROT and CARB, notably:\n"
    ocif+="! - side chain heavy atoms                                                      \n"
    ocif+="                                                                                \n"
    ocif+="define PROT sele ( segid PROA ) end                                             \n"
    ocif+="define CARB sele none end                                                       \n"
    ocif+="define HETE sele ( segid HET* ) end                                             \n"
    ocif+="                                                                                \n"
    ocif+="define BB   sele ( ( type C   .or. type O   .or. -                              \n"
    ocif+="                     type N   .or. type CA  .or. -                              \n"
    ocif+="                     type P   .or. type O1P .or. -                              \n"
    ocif+="                     type O2P .or. type O5' .or. -                              \n"
    ocif+="                     type C5' .or. type C4' .or. -                              \n"
    ocif+="                     type C3' .or. type O3' ) .and. PROT ) .or. -               \n"
    ocif+="                 ( ( type C+ .or. ( type O5 .and. .bonded. type C1 ) .or. -     \n"
    ocif+="                   ( type O6 .and. .bonded. type C2 ) ) .and. CARB ) .or. -     \n"
    ocif+="                 ( .not. hydrogen .and. HETE ) end                              \n"
    ocif+="                                                                                \n"
    ocif+="define SC   sele .not. BB .and. .not. hydrogen .and. ( PROT .or. CARB ) end     \n"
    ocif+="                                                                                \n"
    ocif+="! The following lines are kept for documenting the default restraints           \n"
    ocif+="!cons harm force 1.0 sele BB end                                                \n"
    ocif+="!cons harm force 0.1 sele SC end                                                \n"
    ocif+="                                                                                \n"
    ocif+="! Open output file and write the indexes of the selected residues iteratively   \n"
    ocif+="! The loop structure is based on the CHARMM docs, and uses the selection        \n"
    ocif+="! option .subset.                                                               \n"
    ocif+="                                                                                \n"
    ocif+="open write form name openmm/restraints/prot_pos.txt unit 30                     \n"
    ocif+="echu 30                                                                         \n"
    ocif+="                                                                                \n"
    ocif+="define BB sele BB end                                                           \n"
    ocif+="set bbnum ?nsel                                                                 \n"
    ocif+="set I 0                                                                         \n"
    ocif+="label bblist                                                                    \n"
    ocif+="increment I by 1                                                                \n"
    ocif+="define bbcur sele BB .subset. @I end                                            \n"
    ocif+="calc bbommidx = ?selatom - 1                                                    \n"
    ocif+="if @bbommidx .GE. 0 echo @bbommidx BB ! Avoids -1 BB if sele BB is empty        \n"
    ocif+="if @I .LT. @bbnum goto bblist                                                   \n"
    ocif+="                                                                                \n"
    ocif+="define SC sele SC end                                                           \n"
    ocif+="set scnum ?nsel                                                                 \n"
    ocif+="set I 0                                                                         \n"
    ocif+="label sclist                                                                    \n"
    ocif+="increment I by 1                                                                \n"
    ocif+="define sccur sele SC .subset. @I end                                            \n"
    ocif+="calc scommidx = ?selatom - 1                                                    \n"
    ocif+="if @scommidx .GE. 0 echo @scommidx SC ! Avoids -1 SC if sele SC is empty        \n"
    ocif+="if @I .LT. @scnum goto sclist                                                   \n"
    ocif+="echu                                                                            \n"
    ocif+="                                                                                \n"
    ocif+="! Write out PSF, CRD and PDB                                                    \n"
    ocif+="open write unit 10 card name openmm/step3_input.psf                             \n"
    ocif+="write  psf unit 10 card                                                         \n"
    ocif+="                                                                                \n"
    ocif+="open write unit 10 card name openmm/step3_input.crd                             \n"
    ocif+="write coor unit 10 card                                                         \n"
    ocif+="                                                                                \n"
    ocif+="open write unit 10 card name openmm/step3_input.pdb                             \n"
    ocif+="write coor unit 10 pdb                                                          \n"
    ocif+="                                                                                \n"
    ocif+="stop                                                                            \n"
    ocif+="                                                                                \n"

    ommcinp="${ommcname}.inp"
    if [ -f "${ommcinp}" ]; then
      echo ""
      echo "OpenMM conversion CHARMM input: ${ommcinp} found. Using this."
    else
      printf "%b" "${ocif}" > "${ommcinp}"
      echo ""
      echo "An OpenMM conversion CHARMM input file was not provided."
      echo "Using the standard/built-in scheme. Created file: ${ommcinp}"
      echo "This file can be modified and as long as it is present,"
      echo "it will be used if this process option is selected again."
    fi

    # Create OpenMM coordinates by shifting CHARMM coordinates
    echo ""
    echo "Shifting coordinates to OpenMM origin..."
    echo ""
    charmm_exec < "${ommcinp}" > "${ommcname}.out" &&\
        tail -n 6 ${ommcname}.out &&\
        runmenu

###############################################################################

  elif [ "$menuoption" = "7" ]; then
    echo ""

    # CHARMM-GUI CHARMM Interpreter Scripts for OpenMM
    omm_run="openmm_run.py"
    omm_readinputs="omm_readinputs.py"
    omm_readparams="omm_readparams.py"
    omm_restraints="omm_restraints.py"
    omm_vfswitch="omm_vfswitch.py"
    omm_rewrap="omm_rewrap.py"
    omm_barostat="omm_barostat.py"
    omm_step4_eq="step4_equilibration.inp"
    omm_step5_prod="step5_production.inp"

    cgommpyloc_def=$(realpath "$cgommpyloc_def")

    echo "The default location of the CG OpenMM CHARMM interpreter scripts is:"
    echo "$cgommpyloc_def"
    echo "Press Enter to copy scripts from this location."
    echo "Press d     to provide your own directory."
    echo "Existing scripts in the destination openmm/ will not be overwritten."
    read -r -s -n 1 key

    if [ "$key" = "" ]; then
      cgommpyloc=$cgommpyloc_def
    elif [ "$key" = "d" ]; then
      echo "Enter new directory"
      read cgommpyloc
      if [ "$cgommpyloc" = "" ]; then
        runmenu
      else
        cgommpyloc=$(realpath "$cgommpyloc")
      fi
    else
      runmenu
    fi  

    declare -a scripts=("$omm_run"
                        "$omm_readinputs"
                        "$omm_readparams"
                        "$omm_restraints"
                        "$omm_vfswitch"
                        "$omm_rewrap"
                        "$omm_barostat"
                        "$omm_step4_eq"
                        "$omm_step5_prod")

    #arraylength=${#array[@]}
    #for (( i=0; i<${arraylength}; i++ ));
    echo ""
    for scr in "${scripts[@]}"
    do
      if [ ! -f openmm/"$scr" ]; then
        if [ -f "$cgommpyloc"/"$scr" ]; then
          cp "$cgommpyloc"/"$scr" openmm/ &&\
          echo "Copied "$scr"" 
        else
          echo "File "$cgommpyloc"/"$scr" not found!"
          runmenu
        fi
      fi
    done  
    runmenu
###############################################################################

  elif [ "$menuoption" = "8" ]; then
    echo ""

    # Simple SLURM Script based on CG OpenMM scripts
    suif=''
    suif+='#!/bin/bash                                                                     \n'
    suif+='#SBATCH -p lgpu                                                                 \n'
    suif+='#SBATCH --gres=gpu                                                              \n'
    suif+='#SBATCH --exclude="n[0001-0019]"                                                \n'
    suif+='#SBATCH -J systemlabel                                                          \n'
    suif+='                                                                                \n'
    suif+='source ~/miniconda3/etc/profile.d/conda.sh                                      \n'
    suif+='conda activate openmm                                                           \n'
    suif+='                                                                                \n'
    suif+='init=step3_input                                                                \n'
    suif+='equi_prefix=step4_equilibration                                                 \n'
    suif+='prod_prefix=step5_production                                                    \n'
    suif+='prod_step=step5                                                                 \n'
    suif+='                                                                                \n'
    suif+='# Equilibration                                                                 \n'
    suif+='e_run_param=("-u" "openmm_run.py" "-i" "${equi_prefix}.inp"                     \n'
    suif+='             "-t" "toppar.str" "-p" "${init}.psf"                               \n'
    suif+='             "-c" "${init}.crd" "-b" "sysinfo.dat"                              \n'
    suif+='             "-orst" "${equi_prefix}.rst" "-odcd" "${equi_prefix}.dcd"          \n'
    suif+='            )                                                                   \n'
    suif+='python "${e_run_param[@]}" > ${equi_prefix}.out                                 \n'
    suif+='                                                                                \n'
    suif+='# Production                                                                    \n'
    suif+='p_run_param=("-u" "openmm_run.py" "-i" "${prod_prefix}.inp"                     \n'
    suif+='             "-t" "toppar.str" "-p" "${init}.psf"                               \n'
    suif+='             "-c" "${init}.crd" "-irst" "${equi_prefix}.rst"                    \n'
    suif+='             "-orst" "${prod_prefix}" "-odcd" "${prod_prefix}.dcd"              \n'
    suif+='            )                                                                   \n'
    suif+='python "${p_run_param[@]}" > ${prod_prefix}.out                                 \n'

    submitinp="submit.sh"
    if [ -f "openmm/${submitinp}" ]; then
      echo "SLURM submit script ${submitinp} found."
    else
      printf "%b" "${suif}" > openmm/${submitinp}
      echo "No SLURM submit script provided. Created ${submitinp}"
    fi

    # Edit submit.sh script
    # Update system name - get system name from user
    echo "Input system base name/code"
    read systemname
    sed -e "s/#SBATCH -J .*/#SBATCH -J ${systemname}/g" openmm/submit.sh > openmm/submit_rn.sh && mv openmm/submit_rn.sh openmm/submit.sh
    
    # Update the path to the input PSF 
    sed -e "s/init=.*/init=step3_input/g" openmm/submit.sh > openmm/submit_rn.sh && mv openmm/submit_rn.sh openmm/submit.sh
   
    echo "Enter simulation length in ns (default: 1 ns) and time step in fs (default: 2 fs)"
    echo "Example: 100 2"
    read length dt

    # Evaluate input
    if [[ ${length} == "" && ${dt} == "" ]]; then
        length=500000 # 1 ns
        dt=0.002     # with this timestep
        nstep=`echo "1000000 * ${length} / ${dt}" | bc`
    elif [[ ${length} != "" && ${dt} == "" ]]; then
        dt=0.002
        nstep=`echo "1000000 * ${length} / ${dt}" | bc`
    elif [[ ${length} != "" && ${dt} != "" ]]; then
        nstep=`echo "1000000 * ${length} / ${dt}" | bc`
    fi
    
    echo "Simulation length: ${length} ns"
    echo "Timestep:              ${dt} fs"
    echo "Number of steps:    ${nstep} steps"

   

    # OpenMM Equilibration and Production Input Templates
    ommeqinp="step4_equilibration.inp"
    ommprodinp="step5_production.inp"

    eqif=''
    eqif+='mini_nstep  = 5000                              # Number of steps for minimization                                        \n'
    eqif+='mini_Tol    = 100.0                             # Minimization energy tolerance                                           \n'
    eqif+='                                                                                                                          \n'
    eqif+='gen_vel     = yes                               # Generate initial velocities                                             \n'
    eqif+='gen_temp    = 303.15                            # Temperature for generating initial velocities (K)                       \n'
    eqif+='                                                                                                                          \n'
    eqif+='nstep       = 125000                            # Number of steps to run                                                  \n'
    eqif+='dt          = 0.001                             # Time-step (ps)                                                          \n'
    eqif+='                                                                                                                          \n'
    eqif+='nstout      = 1000                              # Writing output frequency (steps)                                        \n'
    eqif+='nstdcd      = 5000                              # Writing coordinates trajectory frequency (steps)                        \n'
    eqif+='                                                                                                                          \n'
    eqif+='coulomb     = PME                               # Electrostatic cut-off method                                            \n'
    eqif+='ewald_Tol   = 0.0005                            # Ewald error tolerance                                                   \n'
    eqif+='vdw         = Force-switch                      # vdW cut-off method                                                      \n'
    eqif+='r_on        = 1.0                               # Switch-on distance (nm)                                                 \n'
    eqif+='r_off       = 1.2                               # Switch-off distance (nm)                                                \n'
    eqif+='                                                                                                                          \n'
    eqif+='temp        = 303.15                            # Temperature (K)                                                         \n'
    eqif+='fric_coeff  = 1                                 # Friction coefficient for Langevin dynamics                              \n'
    eqif+='                                                                                                                          \n'
    eqif+='pcouple     = no                                # Turn on/off pressure coupling                                           \n'
    eqif+='                                                                                                                          \n'
    eqif+='cons        = HBonds                            # Constraints mehtod                                                      \n'
    eqif+='                                                                                                                          \n'
    eqif+='rest        = yes                               # Turn on/off restraints                                                  \n'
    eqif+='fc_bb       = 400.0                             # Positional restraint force constant for protein backbone (kJ/mol/nm^2)  \n'
    eqif+='fc_sc       = 40.0                              # Positional restraint force constant for protein side-chain (kJ/mol/nm^2)\n'

    prif=''
    prif+='nstep       = 500000                            # Number of steps to run                                                  \n'
    prif+='dt          = 0.002                             # Time-step (ps)                                                          \n'
    prif+='                                                                                                                          \n'
    prif+='nstout      = 1000                              # Writing output frequency (steps)                                        \n'
    prif+='nstdcd      = 50000                             # Writing coordinates trajectory frequency (steps)                        \n'
    prif+='                                                                                                                          \n'
    prif+='coulomb     = PME                               # Electrostatic cut-off method                                            \n'
    prif+='ewald_Tol   = 0.0005                            # Ewald error tolerance                                                   \n'
    prif+='vdw         = Force-switch                      # vdW cut-off method                                                      \n'
    prif+='r_on        = 1.0                               # Switch-on distance (nm)                                                 \n'
    prif+='r_off       = 1.2                               # Switch-off distance (nm)                                                \n'
    prif+='                                                                                                                          \n'
    prif+='temp        = 303.15                            # Temperature (K)                                                         \n'
    prif+='fric_coeff  = 1                                 # Friction coefficient for Langevin dynamics                              \n'
    prif+='                                                                                                                          \n'
    prif+='pcouple     = yes                               # Turn on/off pressure coupling                                           \n'
    prif+='p_ref       = 1.0                               # Pressure (Pref or Pxx, Pyy, Pzz; bar)                                   \n'
    prif+='p_type      = isotropic                         # MonteCarloBarostat type                                                 \n'
    prif+='p_freq      = 100                               # Pressure coupling frequency (steps)                                     \n'
    prif+='                                                                                                                          \n'
    prif+='cons        = HBonds                            # Constraints mehtod                                                      \n'
    prif+='                                                                                                                          \n'
    prif+='rest        = no                                # Turn on/off restraints                                                  \n'

    if [ -f openmm/"${ommeqinp}" ]; then
      echo "OpenMM equilibration input openmm/${ommeqinp} found" 
    else
      printf "%b" "${eqif}" > openmm/${ommeqinp}
      echo "No OpenMM equilibration input provided. Created openmm/${ommeqinp}"
    fi

    if [ -f openmm/"${ommprodinp}" ]; then
      echo "OpenMM production input openmm/${ommprodinp} found" 
    else
      printf "%b" "${prif}" > openmm/${ommprodinp}
      echo "No OpenMM production input provided. Created openmm/${ommprodinp}"
    fi

    # Update the openmm/step5_production.inp with eventual arguments
    # from the user

    # Dirty handling of bc result with many trailing zeros and lacking a prepended zero
    dtps=`echo "${dt} / 1000" | bc -l`
    dtps="0${dtps:0:4}"
    `sed -e "s/nstep       = .*/nstep       = ${nstep}                          # Number of steps to run/g" openmm/step5_production.inp > openmm/step5_production_tmp.inp`
    `sed -e "s/dt          = .*/dt          = ${dtps}                             # Time-step (ps)/g" openmm/step5_production_tmp.inp > openmm/step5_production.inp`
    `/bin/rm openmm/step5_production_tmp.inp`
    
    # Get number of replicas from user
    echo "Enter number of replicas"
    read numberofreplicas

    # Make replica dirs, copy files and edit submit scripts
    curdir=${PWD##*/}
    for rep in $(seq 1 ${numberofreplicas}); do
      mkdir -p ../${curdir}_${rep}
      cp -r openmm/* ../${curdir}_${rep}/
      sed -e "s/#SBATCH -J ${systemname}/#SBATCH -J ${systemname}-${rep}/g"\
      ../${curdir}_${rep}/submit.sh > ../${curdir}_${rep}/submit_rn.sh &&\
      mv ../${curdir}_${rep}/submit_rn.sh ../${curdir}_${rep}/submit.sh &&\
      echo "Created replica at ../${curdir}_${rep}."
    done

    runmenu

###############################################################################

  elif [ "$menuoption" = "a" ]; then
    run_a_submenu

###############################################################################

  elif [ "$menuoption" = "t" ]; then
    run_t_submenu

###############################################################################
  elif [ "$menuoption" = "q" ];then
    echo " <- Goodbye"
    exit 0       

  else
    timeout 2s echo " <- Unrecognized option. Press any key."
    #read -n 1 k <&1
    runmenu
  fi
}

# This builds the main menu and routes the user to the function selected.

runmenu

