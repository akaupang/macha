#!/bin/bash
# For paste toggle in vim
# set pastetoggle=<F2>

runmenu () {
  echo ""
  echo "Manual CHARMM System Setup"
  echo ""
  echo "1: Copy custom system files to current folder (reset)"
  echo "2: Step 1   CHARMM-GUI PDB Reader (modified)"
  echo "3: Step 2.1 CHARMM-GUI Waterbox"
  echo "4: Step 2.2 CHARMM-GUI Ions"
  echo "5: Step 2   CHARMM-GUI Solvator"
  echo "6: Step 3   CHARMM-GUI PBC Setup"
  echo "7: Step 3.1 Convert CHARMM system to an OpenMM system"
  echo "8: Step 3.2 Copy CG OpenMM CHARMM interpreter python scripts to openmm/"
  echo "9: Step 3.3 Set system base name"
  echo "            Set simulation length (and timestep)"
  echo "            Create replicas"
  echo ""
  echo "Press q to quit"
  read -n 1 -p "Process options: " "menuoption"

###############################################################################

  if [ "$menuoption" = "1" ]; then
    echo ""
    echo "Are you sure? (if yes, press 1)"
    read surefire
    if [ "$surefire" == "1" ]; then
      rsync -rtp --exclude '.git' /scratch/data/asmund/repos/mod_systems/* . &&\
      echo "Finished copying" && runmenu
    else
      runmenu
    fi

###############################################################################

  elif [ "$menuoption" = "2" ]; then
    echo ""
    charmm < step1_pdbreader.inp > step1_pdbreader.out && tail -n 6 step1_pdbreader.out && sed -n '/START_PAR/,/END_PAR/p' step1_pdbreader.out > step1_used_parameters.dat
    runmenu
  elif [ "$menuoption" = "3" ]; then
    echo ""
    charmm < step2.1_waterbox.inp > step2.1_waterbox.out && tail -n 6 step2.1_waterbox.out && runmenu

  elif [ "$menuoption" = "4" ]; then
    echo ""
    charmm < step2.2_ions.inp > step2.2_ions.out && tail -n 6 step2.2_ions.out && runmenu

  elif [ "$menuoption" = "5" ]; then
    echo ""
    charmm < step2_solvator.inp > step2_solvator.out && tail -n 6 step2_solvator.out && runmenu

  elif [ "$menuoption" = "6" ]; then
    echo ""
    charmm < step3_pbcsetup.inp > step3_pbcsetup.out && tail -n 6 step3_pbcsetup.out && runmenu

  elif [ "$menuoption" = "7" ]; then
    echo ""

    # Handle an eventual old openmm folder, to allow backup and continued modification (rerun of step3.1). Also handle the case of a new folder (files were recopied, system reset).
    if [ -d "openmm" ]; then
      echo "Found existing OpenMM directory (openmm/)"
      
      # Look for previously backed up folders and give new backup a higher number
      podircount="$(find ./prev_openmm* -maxdepth 0 -mindepth 0 -type d -printf x 2>/dev/null | wc -c)"
      podirnum=`echo "1 + $podircount" | bc`
      podir="prev_openmm_$podirnum"
      mkdir -p "$podir/"

      mv openmm/* "$podir/"
      echo "Backed it up as $podir/"

      # Check if a template directory is present
      if [ -d "new_openmm" ]; then
        `/bin/rm -rf openmm/`
        mv new_openmm/ openmm/
        mkdir -p openmm/restraints/
        echo "Created an OpenMM directory (openmm/) from template files"
      else
        mkdir -p openmm/
        mkdir -p openmm/restraints/
        echo "Created a new OpenMM directory for the current system"
      fi
    else
      if [ -d new_openmm ]; then
        mv new_openmm/ openmm/
        echo "Created an OpenMM directory (openmm/) from template files"
      else
        mkdir -p openmm
        echo "Created a new OpenMM directory for the current system"
      fi
      mkdir -p openmm/restraints/
    fi

    # Transform the CHARMM-GUI toppar.str to an OpenMM toppar.str
    # at least with the formers current formatting
    # The following sed command searches for anything + "toppar/" 
    # and prints out "../" + the rest of the line on which it found 
    # the search parameter, which is now; "line ends with rtf/prm/str".
    `grep -v "!" toppar.str > toppar_tmp.str`

    `sed 's/.* //' toppar_tmp.str | sed -n -e 's/^.*\(rtf$\|prm$\|str$\)/..\/\0/p' > openmm/toppar.str` &&\
    `/bin/rm -f toppar_tmp.str` &&\
    echo "Converted toppar.str (CHARMM format) to openmm/toppar.str (OpenMM format)"


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
    ocif+="echo @bbommidx BB                                                               \n"
    ocif+="if @I .LT. @bbnum goto bblist                                                   \n"
    ocif+="                                                                                \n"
    ocif+="define SC sele SC end                                                           \n"
    ocif+="set scnum ?nsel                                                                 \n"
    ocif+="set I 0                                                                         \n"
    ocif+="label sclist                                                                    \n"
    ocif+="increment I by 1                                                                \n"
    ocif+="define sccur sele SC .subset. @I end                                            \n"
    ocif+="calc scommidx = ?selatom - 1                                                    \n"
    ocif+="echo @scommidx SC                                                               \n"
    ocif+="if @I .LT. @scnum goto sclist                                                   \n"
    ocif+="echu                                                                            \n"
    ocif+="                                                                                \n"
    ocif+="! Write out PSF, CRD and PDB                                                    \n"
    ocif+="open write unit 10 card name openmm/step3.1_omm.psf                             \n"
    ocif+="write  psf unit 10 card                                                         \n"
    ocif+="                                                                                \n"
    ocif+="open write unit 10 card name openmm/step3.1_omm.crd                             \n"
    ocif+="write coor unit 10 card                                                         \n"
    ocif+="                                                                                \n"
    ocif+="open write unit 10 card name openmm/step3.1_omm.pdb                             \n"
    ocif+="write coor unit 10 pdb                                                          \n"
    ocif+="                                                                                \n"
    ocif+="stop                                                                            \n"
    ocif+="                                                                                \n"

    ommcinp="step3.1_omm.inp"
    if [ -f "$ommcinp" ]; then
      echo "OpenMM conversion CHARMM input: $ommcinp found. Using this."
    else
      printf "%b" "$ocif" > "$ommcinp"
      echo "An OpenMM conversion CHARMM input file was not provided."
      echo "Using the standard/built-in scheme. Created file: $ommcinp"
      echo "This file can be modified and as long as it is present,"
      echo "it will be used if this process option is selected again."
    fi

    # Create OpenMM coordinates by shifting CHARMM coordinates
    echo "Shifting coordinates to OpenMM origin..."
    echo ""
    charmm < step3.1_omm.inp > step3.1_omm.out && tail -n 6 step3.1_omm.out && runmenu

###############################################################################

  elif [ "$menuoption" = "8" ]; then
    echo ""

    # CHARMM-GUI CHARMM Interpreter Scripts for OpenMM
    omm_run="openmm_run.py"
    omm_readinputs="omm_readinputs.py"
    omm_readparams="omm_readparams.py"
    omm_restraints="omm_restraints.py"
    omm_vfswitch="omm_vfswitch.py"
    omm_rewrap="omm_rewrap.py"
    omm_barostat="omm_barostat.py"

    cgommpyloc_repo="/scratch/data/asmund/repos/mod_systems/new_openmm"
    cgommpyloc_repo=$(realpath "$cgommpyloc_repo")

    echo "The default location of the CG OpenMM CHARMM interpreter scripts is:"
    echo "$cgommpyloc_repo"
    echo "Press Enter to copy scripts from this location."
    echo "Press d     to provide your own directory."
    echo "Existing scripts in the destination openmm/ will not be overwritten."
    read -r -s -n 1 key

    if [ "$key" = "" ]; then
      cgommpyloc=$cgommpyloc_repo
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
                        "$omm_barostat")

    #arraylength=${#array[@]}
    #for (( i=0; i<${arraylength}; i++ ));
    echo ""
    for scr in "${scripts[@]}"
    do
      if [ ! -f openmm/"$scr" ]; then
        if [ -f "$cgommpyloc"/"$scr" ]; then
          cp "$cgommpyloc"/"$scr" openmm/
          echo "Copied "$scr"" 
        else
          echo "File "$cgommpyloc"/"$scr" not found!"
          runmenu
        fi
      fi
    done  
    runmenu
###############################################################################

  elif [ "$menuoption" = "9" ]; then
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
    if [ -f "openmm/$submitinp" ]; then
      echo "SLURM submit script $submitinp found."
    else
      printf "%b" "$suif" > openmm/$submitinp
      echo "No SLURM submit script provided. Created $submitinp"
    fi

    # Edit submit.sh script
    # Update system name - get system name from user
    echo "Input system base name/code"
    read systemname
    sed -e "s/#SBATCH -J .*/#SBATCH -J ${systemname}/g" openmm/submit.sh > openmm/submit_rn.sh && mv openmm/submit_rn.sh openmm/submit.sh
    
    # Update the path to the input PSF 
    sed -e "s/init=.*/init=step3.1_omm/g" openmm/submit.sh > openmm/submit_rn.sh && mv openmm/submit_rn.sh openmm/submit.sh
   
    echo "Enter simulation length (default: 1 ns) and time step (default: 2 fs)"
    echo "Example: 100 2"
    read nstep dt

    # Evaluate input
    if [[ $nstep == "" && $dt == "" ]]; then
        nstep=500000 # 1 ns
        dt=0.002     # with this timestep
    elif [[ $nstep != "" && $dt == "" ]]; then
        dt=0.002
        nstep=`echo "1000 * $nstep / $dt" | bc`
    elif [[ $nstep != "" && $dt != "" ]]; then
        dt=`echo "$dt/1000" | bc`
        nstep=`echo "1000 * $nstep / $dt" | bc`
    fi
    
    echo "Number of steps: $nstep steps"
    dtfs=`echo "1000 * $dt" | bc`
    echo "Timestep:        $dtfs fs"
   

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

    if [ -f openmm/"$ommeqinp" ]; then
      echo "OpenMM equilibration input openmm/$ommeqinp found" 
    else
      printf "%b" "$eqif" > openmm/$ommeqinp
      echo "No OpenMM equilibration input provided. Created openmm/$ommeqinp"
    fi

    if [ -f openmm/"$ommprodinp" ]; then
      echo "OpenMM production input openmm/$ommprodinp found" 
    else
      printf "%b" "$prif" > openmm/$ommprodinp
      echo "No OpenMM production input provided. Created openmm/$ommprodinp"
    fi

    # Update the openmm/step5_production.inp with eventual arguments
    # from the user
    `sed -e "s/nstep       = .*/nstep       = ${nstep}                            # Number of steps to run/g" openmm/step5_production.inp > openmm/step5_production_tmp.inp`
    `sed -e "s/dt          = .*/dt          = ${dt}                               # Time-step (ps)/g" openmm/step5_production_tmp.inp > openmm/step5_production.inp`
    `/bin/rm openmm/step5_production_tmp.inp`
    
    # Get number of replicas from user
    echo "Enter number of replicas"
    read numberofreplicas

    # Make replica dirs, copy files and edit submit scripts
    curdir=${PWD##*/}
    for rep in $(seq 1 $numberofreplicas); do
      mkdir ../${curdir}_${rep}
      cp -r * ../${curdir}_${rep}/
      sed -e "s/#SBATCH -J ${systemname}/#SBATCH -J ${systemname}-${rep}/g"\
      ../${curdir}_${rep}/openmm/submit.sh > ../${curdir}_${rep}/openmm/submit_rn.sh &&\
      mv ../${curdir}_${rep}/openmm/submit_rn.sh ../${curdir}_${rep}/openmm/submit.sh
    done

    runmenu
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

