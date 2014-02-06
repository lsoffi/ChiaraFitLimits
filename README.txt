channel=category= 0, 1, 2, 3
Gamma_H= 0.1, 3., 6., 10., 15
---------------------- STEP 1: PRODUCE WS AND DATACARD ------------------------
.L ProduceWorkspace.C
ProduceWorkspace(150, 0.1)
ProduceWorkspace(150, 3)
ProduceWorkspace(150, 6)
ProduceWorkspace(150, 10)
ProduceWorkspace(150, 10)

1)This step creates datacard for combine for each channel  and for each
width. They are written in 

----> datacard/HighMass-hgg_8TeV_w*_channel*.txt

2)This step produces also ws used in these datacard and written in:

-----> workspaces/HighMass-hgg.inputbkg.root (BKG, the same for all signals)
-----> worlspace/HighMass-hgg_8TeV_w*.inputsig.root

---------------------- STEP 2: COMBINE DATACARD ------------------------
source combineAllDatacards.csh

This step combine the different channels for each width and wrtie the combined datacard in:

----->datacard/HighMass-hgg_8TeV_w*_COMB.txt

FIXME: togliere ./datacard/ dal path delle shape (va sistemato nello script)

---------------------- STEP 3: RUN COMBINE AND PRODUCE ROOTFILES ------------------------

source combineAllMasses_CH*_w*.csh 
mv Higgs*root CH*/
rm roostat*
This step runs combine for each channel, each width and each mass.

--> CH*/higgsCombine_w*.Asymptotic.mH*.123456.root

source combineAllMasses_COMB_w*.csh
mv Higgs*root COMB/
rm roostat*


---------------------- STEP 4: PLOT LIMITS 1D ------------------------

source makeAllLimitPlots.csh

---------------------- STEP 4: PLOT LIMITS 2D ------------------------
.L 2Dplot.cc+
make2Dplot(0)
make2Dplot(1)
make2Dplot(2)
make2Dplot(3)
make2Dplot(4)<- this make combined limit 2d plot