#!/bin/bash


#
# Simple script to execute ALL commands needed to
# generate the SMAUG isolated galaxy (WLM) ICs
# in the format needed to be read in for Enzo.
#
# This assumes you already have MakeDiskGalaxy compiled
# and loadIC.c compiled. loadIC.c is an interpreter written
# by and obtained from Chia-Yu Hu to parse the MakeDiskGalaxy
# output that has been edited to conform to how Enzo needs the
# data in ProblemType_Agora.
#
# Final data products are zipped together to be moved elsewhere
# if needed (see last command)

MDGNAME="WLM_1E7DM"
MDGINPUT="input_$MGDNAME.txt"
RUN="./"
MDPATH="./../"
ICPATH="./../"

# run make disk galaxy with input
$MDPATH/MakeDiskGalaxy.exe $MDGINPUT > "MakeDiskGalaxy.out"

# extract information from files
$ICPATH/loadIC.exe "ics_$MDGNAME.dat" > "loadIC.out"

# write vcirc
python "$RUN/extract_vcirc.py"

# zip together:
tar -cvf "$MDGNAME_IC.tar" "./disk.dat" "./halo.dat" "./vcirc.dat"
