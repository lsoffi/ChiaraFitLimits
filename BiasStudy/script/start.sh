
#!/bin/bash
#append=local

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

taskName=test


fixed=0
noPois=1
onlyB=0
nToy=10
nJobs=100

muGen=0
jobNumb=0
queue=$1

let nJobs_=$nJobs-1


mkdir -p log/
mkdir -p log/$taskName
mkdir -p batch/
mkdir -p batch/${taskName}
echo "################ Task ${taskName} ##############"    

#DiJet

for jobNumb in `seq 0 $nJobs_`
  do
  for fitFun in Expol Expol2 DiJet
    do
    for genFun in RooKey #Expol Expol2 DiJet  
      do
      for mass in 150 200 250 300 400
	do
        for cat in 0 1 2 3 
	 do
	echo "Processing jobNumb: $jobNumb, fitFun: $fitFun, genFun: $genFun, mass: $mass, cat: $cat" > /dev/stdout    
	
#	analysis_=`echo $PWD | sed 's|/batch||'`
#	analysis=`basename $analysis_`
	
	sample='g'$genFun'_f'$fitFun'_m'$mass'_cat'$cat'_j'$jobNumb'_mu'$muGen
	jobname=$taskName"_"$sample
#	stderr_file=log/$taskname/`basename $sample_ .list`-err.log
	stdout_file=`pwd`/log/$taskName/`basename ${sample}_ .list`-out.log
	#bsub -u $USER@cern.ch -L /bin/bash -J $jobname -q $queue -eo $stderr_file -oo $stdout_file $analysis.sh 
	if [ "$queue" == "local" ]; then  
	    ./script/TestFit.sh $nToy $jobNumb $mass $cat $genFun $fitFun $muGen $fixed $noPois ${onlyB} | tee log/${taskName}/$mass-$genFun-$fitFun-$jobNumb'_'$nJobs-$muGen.log
        exit 0
	else
           
	    cat <<EOF >  batch/$taskName/${jobname}.sh
#!/bin/sh
`pwd`/script/TestFit.sh $nToy $jobNumb $mass $cat $genFun $fitFun $muGen $fixed $noPois ${onlyB}
EOF
    echo "bsub -J $jobname -q $queue -o $stdout_file < `pwd`/batch/$taskName/${jobname}.sh"
	    bsub -J $jobname -q $queue -o $stdout_file < `pwd`/batch/$taskName/${jobname}.sh
	fi
      done
    done
  done
done
done
	
