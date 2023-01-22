[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/akaupang/macha/workflows/CI/badge.svg)](https://github.com/akaupang/macha/actions)
[![codecov](https://codecov.io/gh/wiederm/transformato/branch/master/graph/badge.svg)](https://codecov.io/gh/akaupang/macha/branch/main)


## maNUAL chaRMM (macha)

#### Installation

macha is available as a conda package and can be installed in any conda environment with:

`conda install -c johanneskarwou -c conda-forge macha`


#### Note:
Please note that this an early incarnation, which may or may not be suitable for general use.
It is recommended to set an alias for the script for convenient usage, e.g. in your .bashrc:
\# Make sure that the charmm binary is in the $PATH or you have set a location in the script.
export PATH="/path/to/charmm/bin:$PATH"
\# Create an alias for the script e.g. in your .bashrc .
alias macha="/path/to/manual_charmm_system_setup.sh"
