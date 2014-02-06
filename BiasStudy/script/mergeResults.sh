#!/bin/bash


mass=300
#sigGen=0
cat=0
perc=0
winSize=10

  for fitFun in  ExpPAR  # DiJetEXP DiJetEXPOL ExpolPL
    do
    for genFun in Expol DiJet #
      do
      for mass in  150 #250 350 450 550 650 750 850
	do
	for cat in  0 #1 2 3 #1 2 3 
	 do
	 for perc in 10 #0  #5 10  #10 #5 10
	  do
          for suffix in FIX
            do
	   echo   "hadd -f TestSWResults/mergeDir/biasCheck-sigGen_0-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_${suffix}.root TestSWResults/biasCheck-sigGen_0-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*_${suffix}.root"
	      
	      hadd -f TestSWResults/mergeDir/biasCheck-sigGen_1-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_${suffix}.root TestSWResults/biasCheck-sigGen_842-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*_${suffix}.root
	     # cp TestSWResults/mergeDir/biasCheck-sigGen_0-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root /afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/
	      echo "Done mass: $mass cat: $cat gen: $genFun fit: $fitFun  %: $perc range: $suffix" 
   
	   
	 done
	done 	
      done  
    done
done 
done


#mass450
#250_1000 250_800 250_700 250_600 300_1000 300_800 300_700 300_600 130_1000 130_800 130_700 130_600 200_1000 200_800 200_700 200_600

#mass250
#130_350 130_450 130_550 150_350 150_450 150_550 200_350 200_450 200_550

#650
#350_1000 350_900 350_800 450_1000 450_900 450_800 550_1000 550_900 550_800

#m150
#130_200 130_230 130_250

#350
#200_450 200_550 200_650 250_450 250_550 250_650 300_450 300_550 300_650