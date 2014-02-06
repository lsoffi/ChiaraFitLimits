#!/bin/bash



case $# in 
    11)
	nToy=$1
	jobNumb=$2
	mass=$3
	cat=$4
	win=$5
	genFun=$6
	fitFun=$7
	sGen=$8
	perc=$9	
	onlyB=$10
	suf=${11}
       	taskName=${12}
	;;
    *)
	echo "Wrong number of parameters"
	exit 1
	;;
esac
fixed=0
nJobs=1
plotStep=1
noPois=0

echo $suf
echo $perc
sample='g'$genFun'_f'$fitFun'_m'$mass'_w'$perc'_cat'$cat'_nSig'$sGen'_j'$jobNumb'_'$suf
taskName=$sample
echo " ***** ${taskName} ********* "
#export local_output_dir=./output
export xrootd_server=pccmsrm27
export output_dir=/cms/local/soffi/BiasStudy

#CMSSW_DIR=/afs/cern.ch/cms/CAF/CMSPHYS/PHYS_EGAMMA/electrons/meridian/CMSSW425/src
CMSSW_DIR=/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy
BASE_DIR=/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy

#shapeFile=${BASE_DIR}/events_45.0_25.0_-10000.0_-10000.0_-10000.0_-10000.0_70.0_0.0_0.0_-10000.0_0.0_0.0_0.0_0.0_0.0_300.0_0.0_-1_0_0_1_4_cs_metTagWS.root
#shapeFile=${BASE_DIR}/events_45.0_25.0_-10000.0_-10000.0_-10000.0_-10000.0_60.0_0.0_0.0_-10000.0_0.0_0.0_0.0_0.0_0.0_300.0_0.0_-1_0_0_1_4_cs_w_metTagWS.root
#shapeFile=${BASE_DIR}/events_met70_2012_ebeb-ws.root
#shapeFile=${BASE_DIR}/HGG_cat7_cs_btagloose.root
#signalShapeFile=${BASE_DIR}/histogg_cat7.root

cd ${CMSSW_DIR}/src
eval `scramv1 runtime -sh`

if [ -n "${WORKDIR:-x}" ]; then
    WORKDIR=${TMPDIR}
fi
cd ${WORKDIR}
echo "Running in directory `pwd` on `uname -a`"


xrd ${xrootd_server} mkdir  ${output_dir}/


mkdir ${sample}
cd ${sample}
#export sample=`dirname $1 | sed 's|.*filelist/||'`
#export index=`basename $1 .list | sed 's|.*-||'`
#dir=`pwd`/batch
#mkdir -p batch
#source $dir/base_script

 mkdir TestSWResults
 mkdir SWplots
 rm SWplots/*png
 rm SWplots/toys/*png
 mkdir SWplots/toys/
 mkdir log
#shapeFile=events_45.0_25.0_-10000.0_-10000.0_-10000.0_-10000.0_60.0_0.0_0.0_-10000.0_0.0_0.0_0.0_0.0_0.0_300.0_0.0_-1_0_0_1_4_cs_metTagWS.root

#shapeFile=events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.5_350.0_2.6_-1_-1_0_4_cs-ws.root
#shapeFile=events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.5_350.0_0.0_-1_-1_0_4_cs-ws-test.root

#options="-t $nToy -p $plotStep -nJ $nJobs -j $jobNumb -m $mass -w $winSize -wide -pG -gen $genFun -fit $fitFun -f $shapeFile -dataf ${signalShapeFile}" 
InputFit=${BASE_DIR}workspaces/HighMass-hgg.m${mass}.0_w${perc}_inputBias_truth_8TeV_$suf.root
InputGen=${BASE_DIR}/workspaces/HighMass-hgg.m${mass}.0_w${perc}_inputBias_truth_8TeV_$suf.root
options="-t $nToy -p $plotStep -nJ $nJobs -j $jobNumb -mass $mass -width $perc -cat $cat -pG -gen $genFun -fit $fitFun  -f_fit $InputFit -f_gen $InputGen -sigGen -nGen $sGen -w $win -su $suf" 

#only B fit
if [ "$onlyB" != "0" ]; then
    onlyBOpt= "" #" -b"
    options=${options}${onlyBOpt}
fi

command="${BASE_DIR}/build/TestSW_wind.exe"


if [ "$fixed" != "0" ]; then
    fixString="" # -fixed"
    options=$options$fixString
fi

echo "$command $options"


#mkdir -p log/${taskName}

$command $options > log/${sample}.log

#mkdir -p output/
#mv TestSWResults output/
#mv SWplots output/
#mv log output/


echo "$sample done"

cd SWplots
tar cvzf ../$sample.tar.gz *png toys/*.png
cd ..

for file in TestSWResults/*.root; do 
    xrd ${xrootd_server} rm ${output_dir}/${sample}/`basename $file`
    #xrdcp $file root://${xrootd_server}/${output_dir}/`basename $file`; 
    cp $file /afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/`basename $file`; 
done

echo "$sample copy done"


exit 0

