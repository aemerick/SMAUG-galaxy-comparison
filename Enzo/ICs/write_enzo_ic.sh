#!/bin/bash

MDGNAME="WLM_1E7DM"
MDGINPUT="input_$MGDNAME.txt"
RUN="/home/aemerick/work/SMAUG_IC"

# run make disk galaxy with input
./../MakeDiskGalaxy.exe $MDGINPUT > "MakeDiskGalaxy.out"

# extract information from files
./../loadIC.exe "ics_$MDGNAME.dat" > "loadIC.out"

# write vcirc
python "$RUN/extract_vcirc.py"

# zip together:
tar -cvf "$MDGNAME_IC.tar" "./disk.dat" "./halo.dat" "./vcirc.dat"
