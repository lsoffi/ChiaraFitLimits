#!/bin/bash

nToy=20
plotStep=1
mass=300
nJobs=1
let nJobs_=$nJobs-1
jobNumb=0
muGen=0
cat=0
perc=$1

mergeDir=TestSWResults/mergeDir
if [ ! -e "$mergeDir" ]; then mkdir $mergeDir; fi



for jobNumb in `seq 0 $nJobs_`
  do
#exp3_3 pow2_2 exp2_2  exp3_2
  for fitFun in  ExpPAR
    do
    for genFun in  ExpPAR
      do
      for mass in 250 #200 250 300 400 
	do
	for cat in 0 #1 2 3 
	 do
         for suffix in 200_500
          do
	  (
	      InputPrefix=/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy
	      InputFit=$InputPrefix/workspaces/HighMass-hgg.m${mass}.0_w${perc}_inputBias_truth_8TeV_${suffix}.root
	      InputGen=$InputPrefix/workspaces/HighMass-hgg.m${mass}.0_w${perc}_inputBias_truth_8TeV_${suffix}.root
	      
	      
	      $InputPrefix/build/TestSW_wind.exe -t $nToy -p $plotStep -nJ $nJobs -j $jobNumb -mass $mass -width $perc -cat $cat -pG \
		  -gen $genFun -fit $fitFun -p $step   \
		  -f_fit $InputFit -f_gen $InputGen \
		   -sigGen -nGen 1  -w 10 -su $suffix\
	         # &> $InputPrefix/log/$mass-$cat-$genFun-$fitFun-$jobNumb'_'$nJobs-$muGen.log
	    echo "Done mass: $mass cat: $cat gen: $genFun fit: $fitFun $jobNumb/$nJobs"
	    
# 		else
# 		    echo "Already done: $mass $genFun $fitFun $jobNumb/$nJobs"
# 		fi
	  ) &
	  done
	done 		# fit fun
#	  wait
      done   # mass
      wait
    done # genFun
    wait
done # jobNumb
done #suffix


