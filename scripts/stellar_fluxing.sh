#!/bin/bash
#script to run stellar_fluxing.py

#output file name
fname=$1
fname="$fname.txt"

#1st argument is masknumber, 2nd argument is boolean whether bluemask or not
python stellar_fluxing.py $1 $2 > ../results/fluxing/delMag/$fname
echo 'Done'