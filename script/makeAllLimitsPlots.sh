#!/bin/bash

for cat in COMB
  do
  for width in  10.00
    do
    echo  $width
   
    python limitPlotter_livia.py -M Asymptotic  -e -v -p /afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/${cat} -c ${cat} -w _w${width} -t _w${width}

  done #{width}
done #mass


