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
nToy=20
nJobs=50

muGen=0

jobNumb=0
queue=$1
perc=0
sigGen=1
let nJobs_=$nJobs-1


mkdir -p log/
mkdir -p log/$taskName
mkdir -p batch/
mkdir -p batch/${taskName}
echo "################ Task ${taskName} ##############"    

#DiJet

for jobNumb in `seq 0 $nJobs_`
  do
  for fitFun in ExpPAR  #DiJetEXP DiJetEXPOL DiJetPL #DiJetEXP ExpolPL 
    do
    for genFun in Expol DiJet #PL DiJetPL DiJetEXPOL DiJetEXP RooKey3 RooKey4
      do
      for mass in 150 250 350 450 550 650 750 850 
	do
        for cat in 0 #1 2 3 
	 do
	  for win in  10
	    do
	    for perc in 10  #5 10  #10 
	      do
	      for suffix in FIX
		do
		
		echo "Processing jobNumb: $jobNumb, fitFun: $fitFun, genFun: $genFun, mass: $mass, perc: $perc cat: $cat, winSize: $win range: $suffix" > /dev/stdout    
		
#	analysis_=`echo $PWD | sed 's|/batch||'`
#	analysis=`basename $analysis_`
		
		sample='g'$genFun'_f'$fitFun'_m'$mass'_w'$perc'_cat'$cat'_winSize'$win'_j'$jobNumb'_nSig'$sigGen_'_'${suffix}
		jobname=$taskName"_"$sample
#	stderr_file=log/$taskname/`basename $sample_ .list`-err.log
		stdout_file=`pwd`/log/$taskName/`basename ${sample}_ .list`-out.log
	#bsub -u $USER@cern.ch -L /bin/bash -J $jobname -q $queue -eo $stderr_file -oo $stdout_file $analysis.sh 
		if [ "$queue" == "local" ]; then  
		    ./script/TestFit_win.sh $nToy $jobNumb $mass $cat $win $genFun $fitFun $sigGen $perc  ${onlyB} ${suffix}| tee log/${taskName}/$mass-$perc-$cat-$win-$genFun-$fitFun-$jobNumb'_'$nJobs-$sigGen'_'${suffix}.log
		    exit 0
		else
		    
		    cat <<EOF >  batch/$taskName/${jobname}.sh
#!/bin/sh
`pwd`/script/TestFit_win.sh $nToy $jobNumb $mass $cat $win $genFun $fitFun $sigGen $perc  ${onlyB} ${suffix}
EOF
		    echo "bsub -J $jobname -q $queue -o $stdout_file < `pwd`/batch/$taskName/${jobname}.sh"
		 bsub -J $jobname -q $queue -o $stdout_file < `pwd`/batch/$taskName/${jobname}.sh
		fi
	      done
	    done
	 done
	done
      done
    done
  done
done	


#m450
#250_1000 250_800 250_700 250_600 300_1000 300_800 300_700 300_600 130_1000 130_800 130_700 130_600 200_1000 200_800 200_700 200_600 

#m250
#130_350 130_450 130_550 150_350 150_450 150_550 200_350 200_450 200_550


#m650
#350_1000 350_900 350_800 450_1000 450_900 450_800 550_1000 550_900 550_800


#m850
#450_1000 550_1000 650_1000


#m150
#130_300 130_250 130_230 130_200


#350
#200_450 200_550 200_650 250_450 250_550 250_650 300_450 300_550 300_650
