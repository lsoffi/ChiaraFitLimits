#!/bin/bash


case $# in 
    6)
	nToy=$1
	jobNumb=$2
	mass=$3
	genFun=$4
	fitFun=$5
	muGen=$6
	;;
    *)
	exit 1
	;;
esac

nJobs=1
plotStep=500
winSize=10

export local_output_dir=./output
export output_dir=$CASTOR_HOME/
#export sample=`dirname $1 | sed 's|.*filelist/||'`
#export index=`basename $1 .list | sed 's|.*-||'`
dir=/afs/cern.ch/user/s/shervin/Higgs/TestFit/batch
#source $dir/base_script

#mkdir TestSWResults
#mkdir SWplots
#mkdir log
#shapeFile=events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_0.0_350.0_0.0_-1_-1_0_4_cs-ws.root
shapeFile=events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.5_350.0_2.6_-1_-1_0_4_cs-ws.root

case $muGen in
    0)
	./build/TestSW.exe -t $nToy -p $plotStep -nJ $nJobs -j $jobNumb -m $mass -w $winSize -wide -pG \
	    -gen $genFun -fit $fitFun  \
	    -f $shapeFile \
	    -dataf histo_massgg_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.5_350.0_2.6_-1_-1_0_4.root &> log/$mass-$genFun-$fitFun-$jobNumb'_'$nJobs-$muGen.log
	    echo "Done $mass $genFun $fitFun $jobNumb/$nJobs"
	    ;;
    *)
	./build/TestSW.exe -t $nToy -p $plotStep -nJ $nJobs -j $jobNumb -m $mass -w $winSize -wide -pG \
	    -gen $genFun -fit $fitFun -sigGen -muGen $muGen  \
	    -f $shapeFile \
	    -dataf histo_massgg_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.5_350.0_2.6_-1_-1_0_4.root &> log/$mass-$genFun-$fitFun-$jobNumb'_'$nJobs-$muGen.log
	;;
esac

#mkdir output/
#mv TestSWResults output/
#mv SWplots output/
#mv log output/

sample=g$genfun'_f'$fitFun'_m'$mass'_j'$jobNumb
echo "$sample done"
#tar -xzf $sample.tar.gz output/
#xrdcp $sample.tar.gz root://eoscms//eos/cms/store/caf/user/shervin/
#rfcp $sample.tar.gz /castor/cern.ch/user/s/shervin/

exit 0

