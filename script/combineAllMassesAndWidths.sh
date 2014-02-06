#!/bin/bash

for mass in 150 #800 #200 300 400 500 600 700 #150 250 350 450 550  650 750 850
  do
  for width in 0.10 10.00
    do
    echo $mass $width
   
    combine  /afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/datacardWithSyst/HighMass-hgg_8TeV_m${mass}.00_w${width}_COMB.txt -M Asymptotic -m ${mass} -D data_obs --run=expected  -U --noFitAsimov --saveToys -S 0 -n _w${width}

  done #{width}
done #mass


