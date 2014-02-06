#!/bin/bash

for mass in 800.00 # 200.00 300.00 400.00 500.00 600.00 700.00 #150.00 250.00 350.00 450.00 550.00  650.00 750.00 850.00
  do
  for width in 10.00 0.10
    do
    echo $mass $width
    rm HighMass-hgg_8TeV_m${mass}_w${width}_COMB.txt
	       
    ../combineCards.py Name0=HighMass-hgg_8TeV_m${mass}_w${width}_channel0.txt Name1=HighMass-hgg_8TeV_m${mass}_w${width}_channel1.txt Name2=HighMass-hgg_8TeV_m${mass}_w${width}_channel2.txt Name3=HighMass-hgg_8TeV_m${mass}_w${width}_channel3.txt > HighMass-hgg_8TeV_m${mass}_w${width}_COMB.txt

  done #{width}
done #mass


