#!/bin/bash

for mass in 150 # 250 350 450 550 650 750 850 
  do
  for width in 10.00 
    do
    echo $mass $width
   
    echo "combine  /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/datacardWithSyst/HighMass-hgg_8TeV_m${mass}.00_w${width}_COMB.txt -M Asymptotic -m ${mass} -D data_obs --run=expected  -U --noFitAsimov --saveToys -S 0 -n _w${width}"

  done #{width}
done #mass


