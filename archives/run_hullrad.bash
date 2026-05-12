#!/bin/bash
#####################################
# Script for running HullRad on a collection of PDB files on UNIX or OSX
#
# Edit the basename of the PDB files below, i.e. "*.pdb" to your basename
#
# Edit the full path of your python installation that has the following libraries:
#       numpy, scipy
# if you are using these instead of the qhull version of qconvex
#####################################

if  [ -e 'Asphericity.dat' ]
then
  rm Asphericity.dat
fi
touch Asphericity.dat

if  [ -e 'Dmax.dat' ]
then
  rm Dmax.dat
fi
touch Dmax.dat

if  [ -e 'Rg.dat' ]
then
  rm Rg.dat
fi
touch Rg.dat

if  [ -e 's.dat' ]
then
  rm s.dat
fi
touch s.dat

if  [ -e 'filenames.lst' ]
then
  rm filenames.lst	
fi
touch filenames.lst

for file in $(ls s*.pdb)
do
  echo $file 
  echo $file >> filenames.lst
  /opt/local/Library/Frameworks/Python.framework/Versions/3.9/bin/python3.9 \
  HullRadV10.py $file > HRout.tmp
  grep 'Rg\(Anhydrous\)' HRout.tmp | awk '{print $3}' >> Rg.dat
  grep 'Dmax ' HRout.tmp | awk '{print $3}' >> Dmax.dat
  grep 'Asphericity' HRout.tmp | awk '{print $3}' >> Asphericity.dat
  grep 's20,w      ' HRout.tmp | awk '{print $3}' >> s.dat
done

rm HRout.tmp
