#!/bin/bash


mass=300
#sigGen=0
cat=0
perc=0
winSize=10

  for fitFun in  DiJetEXPOL # DiJetEXP DiJetEXPOL ExpolPL
    do
    for genFun in  RooKey4 #RooKey4   #RooKey DiJetPL DiJetEXP DiJetEXPOL ExpolPL
      do
      for mass in 150 200 250 300 350 400 450 500 550 650 750 850
	do
	for cat in 0 #1 2 3 
	 do
	 for perc in 0 5 10 #5 10
	  do
	  if [ $mass -eq 150 ]; then
	      echo   "hadd -f TestSWResults/mergeDir/biasCheck-sigGen_57.274-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_57.2743-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root"
	      
	      hadd -f TestSWResults/mergeDir/biasCheck-sigGen_57.274-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_57.2743-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root
	     # cp TestSWResults/mergeDir/biasCheck-sigGen_57.2743-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/
	      echo "Done mass: $mass cat: $cat gen: $genFun fit: $fitFun  %: $perc" 
	  fi	  
	  if [ $mass -eq 200 ]; then
	      echo   "hadd -f TestSWResults/mergeDir/biasCheck-sigGen_73.980-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_73.98-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root"
	      
	      hadd -f TestSWResults/mergeDir/biasCheck-sigGen_73.980-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_73.98-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root
	      #cp TestSWResults/mergeDir/biasCheck-sigGen_73.98-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/
	   echo "Done mass: $mass cat: $cat gen: $genFun fit: $fitFun  %: $perc" 
	  fi
	  if [ $mass -eq 250 ]; then
	      echo   "hadd -f TestSWResults/mergeDir/biasCheck-sigGen_86.126-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_86.1266-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root"
	      
	      hadd -f TestSWResults/mergeDir/biasCheck-sigGen_86.126-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_86.1266-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root
	      #cp TestSWResults/mergeDir/biasCheck-sigGen_86.1266-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/
	      echo "Done mass: $mass cat: $cat gen: $genFun fit: $fitFun  %: $perc"
	  fi
	  if [ $mass -eq 300 ]; then
	      echo   "hadd -f TestSWResults/mergeDir/biasCheck-sigGen_99.099-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_99.099-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root"
	     
	      hadd -f TestSWResults/mergeDir/biasCheck-sigGen_99.099-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_99.099-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root
	     #cp TestSWResults/mergeDir/biasCheck-sigGen_99.099-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/
	     echo "Done mass: $mass cat: $cat gen: $genFun fit: $fitFun  %: $perc"
	  fi
	  if [ $mass -eq 350 ]; then
	      echo   "hadd -f TestSWResults/mergeDir/biasCheck-sigGen_107.403-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_107.403-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root"
	     
	      hadd -f TestSWResults/mergeDir/biasCheck-sigGen_107.403-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_107.403-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root
	      #cp TestSWResults/mergeDir/biasCheck-sigGen_107.403-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/
	      echo "Done mass: $mass cat: $cat gen: $genFun fit: $fitFun  %: $perc"
	  fi
	  if [ $mass -gt 350 ]; then
	      echo   "hadd -f TestSWResults/mergeDir/biasCheck-sigGen_115.707-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_115.707-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root"
	      
	      hadd -f TestSWResults/mergeDir/biasCheck-sigGen_115.707-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_115.707-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root
	     # cp TestSWResults/mergeDir/biasCheck-sigGen_115.707-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/
	      echo "Done mass: $mass cat: $cat gen: $genFun fit: $fitFun  %: $perc"
	  fi
	  
	 #hadd -f TestSWResults/mergeDir/biasCheck-sigGen_0-gDiJetEXPOL_fDiJetPL_m850_w0_cat0_w10.root TestSWResults/biasCheck-sigGen_0-gDiJetEXPOL_fDiJetPL_m850_w0_cat0_w10_j*.root

#	 echo   "hadd -f TestSWResults/mergeDir/biasCheck-sigGen_${sigGen}-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_${sigGen}-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root"

	# hadd -f TestSWResults/mergeDir/biasCheck-sigGen_${sigGen}-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root TestSWResults/biasCheck-sigGen_${sigGen}-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}_j*.root
	# cp TestSWResults/mergeDir/biasCheck-sigGen_${sigGen}-g${genFun}_f${fitFun}_m${mass}_w${perc}_cat${cat}_w${winSize}.root /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/
	#   echo "Done mass: $mass cat: $cat gen: $genFun fit: $fitFun  %: $perc"
	    
	   
	 done
	done 	
      done  
    done
done 



