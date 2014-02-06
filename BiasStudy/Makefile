#Shervin please stop use svn

OPT=-Wall
ROOT_LIB=`root-config --libs --glibs`
ROOT_FLAGS=`root-config --cflags`
ROOT_INCLUDE=`root-config --incdir`


SCRAMTOOL=$(shell scram || echo 1)
ifeq ($(SCRAMTOOL),1)
	ROOFIT_LIB="-lRooFit"
	ROOFIT_INCLUDE="./"
else
	ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
	ROOFIT_LIB := -l$(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIB= | sed 's|LIB=||')
	ROOFIT_LIB += -l$(shell cd $(CMSSW_BASE); scram tool info roofit | grep LIB= | sed 's|LIB=||')
	ROOFIT_LIBDIR = -L$(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')
	ROOFIT_LIB+=$(ROOFIT_LIBDIR)
endif

#	cd ~/CMSSW_3_6_2/src; 
#INCLUDE_ROOFIT=`scram tool info roofitcore | grep INCLUDE | sed 's|INCLUDE=||'`
#INCLUDE_ROOFIT=/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms6/include
#LD_ROOFIT_PATH=/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms6/lib
#LD_ROOFIT_PATH=`scram tool info roofitcore| grep LIBDIR | sed 's|LIBDIR=||'`
#INCLUDE_MLFIT=/afs/cern.ch/user/s/shervin/CMSSW_3_5_8/src/MLFit/workdir/MLFit
#export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$INCLUDE_ROOFIT #:$INCLUDE_MLFIT
#    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LD_ROOFIT_PATH
#alias g++="g++ -Wall `root-config --cflags --libs --glibs` -L$LD_ROOFIT_PATH -lRooFit"

# -L/usr/lib
MYTOOL_DIR=../../myTools
MYTOOL_OBJ=$(MYTOOL_DIR)/*/obj/*.o
MYTOOL_INC=$(MYTOOL_DIR)/*/include/

INCLUDEDIR=./include
SRCDIR=./src
BUILDDIR=./build
OBJDIR=./lib

OBJ_DIR=./obj
BUILD=./build

#MODULES=hist_class
#LINKER=g++
INCLUDE=-I$(INCLUDEDIR) -I$(ROOT_INCLUDE)  -I$(ROOFIT_INCLUDE)
#-I$(MYTOOL_DIR)/trigger_class/include/ -I$(MYTOOL_DIR)/json_class/include

default: test 

test: $(BUILD)/TestSW.exe $(BUILD)/TestSW_wind.exe $(BUILD)/TestSW_noPois.exe $(BUILD)/TestSW_wind_noPois.exe 
$(BUILD)/TestSW.exe: $(SRCDIR)/TestSW.cc 
	g++ $(OPT) $(INCLUDE) $(ROOT_LIB) $(ROOFIT_LIB) -o $(BUILD)/TestSW.exe $(SRCDIR)/TestSW.cc
$(BUILD)/TestSW_wind.exe: $(SRCDIR)/TestSW_wind.cc 
	g++ $(OPT) $(INCLUDE) $(ROOT_LIB) $(ROOFIT_LIB) -o $(BUILD)/TestSW_wind.exe $(SRCDIR)/TestSW_wind.cc
$(BUILD)/TestSW_noPois.exe: $(SRCDIR)/TestSW.cc 
	g++ $(OPT) -DNOPOIS $(INCLUDE) $(ROOT_LIB) $(ROOFIT_LIB) -o $(BUILD)/TestSW_noPois.exe $(SRCDIR)/TestSW.cc
$(BUILD)/TestSW_wind_noPois.exe: $(SRCDIR)/TestSW_wind.cc 
	g++ $(OPT) -DNOPOIS $(INCLUDE) $(ROOT_LIB) $(ROOFIT_LIB) -o $(BUILD)/TestSW_wind_noPois.exe $(SRCDIR)/TestSW_wind.cc



#$(MYTOOL_DIR)/functions/obj/setTDRStyle.o








mytool: 
	cd $(MYTOOL_DIR); make;	cd -

