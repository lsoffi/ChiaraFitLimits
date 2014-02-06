
#!/bin/bash
#append=local
queue=8nm

usage(){
    echo "Usage: `basename $0` joblist_file" > /dev/stderr
}

case $# in 
    0)
	echo -n "Error: " > /dev/stderr
	usage
	exit 1
	;;
    *)
	;;
esac
muGen=$1
nToy=1000
nJobs=1
let nJobs_=$nJobs-1
jobNumb=$2
#for jobNumb in `seq 0 $nJobs_`
#  do
#exp3_3 pow2_2 exp2_2  exp3_2
#1exp 1pow 2pol pow2_1 exp2_1 exp3_1 1pol  3pol 4pol 
  for fitFun in  3pol
    do
    for genFun in 1pow 2pol 1exp
      do
      for mass in 110 115 120 130 140
	do

	echo "Processing jobNumb: $jobNumb, fitFun: $fitFun, genFun: $genFun, mass: $mass" > /dev/stdout    
	
	analysis_=`echo $PWD | sed 's|/batch||'`
	analysis=`basename $analysis_`
	
	sample=g$genfun'_f'$fitFun'_m'$mass'_j'$jobNumb'_mu'$muGen
	jobname=$analysis"_"$sample
	stderr_file=log/$sample/$queue/`basename $sample_ .list`-err.log
	stdout_file=log/$sample/$queue/`basename $sample_ .list`-out.log
#  echo $i
#  echo $sample
      #jobname=
#      echo "bsub -u meridian@cern.ch -L /bin/bash -J $jobname -q $queue -eo $stderr_file -oo $stdout_file $analysis.sh $i"
	#bsub -u $USER@cern.ch -L /bin/bash -J $jobname -q $queue -eo $stderr_file -oo $stdout_file $analysis.sh 
	./script/TestFit.sh $nToy $jobNumb $mass $genFun $fitFun $muGen
	done
      done
  done
#done

	
