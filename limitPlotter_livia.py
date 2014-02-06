#!/usr/bin/env python
# Original Authors - Nicholas Wardle, Nancy Marinelli, Doug Berry

# Major cleanup from limit-plotter-complete.py
#-------------------------------------------------------------------------
# UserInput
from optparse import OptionParser
parser=OptionParser()
parser.add_option("-c", "--category", dest="cat", default="0", type="str")
parser.add_option("-t", "--suffix", dest="suffix", default="0", type="str")
parser.add_option("-w","--width", dest="width", default="0.10", type="str")
parser.add_option("-M","--Method",dest="Method")
parser.add_option("-r","--doRatio",action="store_true")
parser.add_option("-s","--doSmooth",action="store_true")
parser.add_option("-b","--bayes",dest="bayes")
parser.add_option("-o","--outputLimits",dest="outputLimits")
parser.add_option("-e","--expectedOnly",action="store_true")
parser.add_option("-p","--path",dest="path",default="",type="str")
parser.add_option("-v","--verbose",dest="verbose",action="store_true")
parser.add_option("-a","--append",dest="append",default="",help="Append string to filename")
parser.add_option("","--sideband",dest="sideband",default=False,action="store_true")
parser.add_option("","--addline",action="append",type="str",help="add lines to the plot file.root:color:linestyle:legend entry", default = [])
parser.add_option("","--show",action="store_true")
parser.add_option("","--pval",action="store_true")
parser.add_option("","--addtxt",action="append",type="str", help="Add lines of text under CMS Preliminary",default=[])
parser.add_option("","--square",dest="square",help="Make square plots",action="store_true")
parser.add_option("","--nogrid",dest="nogrid",help="Remove grid from plots",action="store_true")

(options,args)=parser.parse_args()

# Standard Imports and calculators
import ROOT
import array,sys,numpy
ROOT.gROOT.ProcessLine(".x /afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/tdrstyle.cc")

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
www_path="/afs/cern.ch/user/s/soffi/www/Limits"
#---------open out file txt-----
out_file=open("limitsResults_"+options.cat+options.suffix+".txt", "a")
#-------------------------------------------------------------------------
# Configuration for the Plotter
OBSmasses = []
EXPmasses = []

if options.width == "_w0.10":OBSmassesT = [ 150,160, 170, 180, 190, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 500, 530,560,  630,  675, 695,700, 730,760,800, 830, 850]
if options.width == "_w2.00":OBSmassesT = [ 150,160, 170, 180, 190, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 510,520, 530, 560, 590, 640, 670,  710,  730,760,800,  850]
if options.width == "_w5.00":OBSmassesT = [ 150,160, 170, 180, 190, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 500, 530,560,570, 580, 590, 620, 650, 660,700, 730,760,800,  850]
if options.width == "_w7.00":OBSmassesT = [ 150,160, 170, 180, 190, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 500, 530,570, 580,610,  630, 660,710, 750, 780,800, 830, 850]
if options.width == "_w10.00":OBSmassesT = [ 150,160, 170, 180, 190, 200,220,  260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 500, 530,560, 600,610,640, 660, 750, 770,800, 850]

if options.width == "_w0.10":EXPmassesT = [ 150,160, 170, 180, 190, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 500, 530,560,  630,  675, 695,700, 730,760,800, 830, 850]
if options.width == "_w2.00":EXPmassesT = [ 150,160, 170, 180, 190, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 510,520, 530, 560, 590, 640, 670, 710,  730,760,800,  850]
if options.width == "_w5.00":EXPmassesT = [ 150,160, 170, 180, 190, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 500, 530,560,570, 580, 590, 620, 650, 660,700, 730,760,800,  850]
if options.width == "_w7.00":EXPmassesT = [ 150,160, 170, 180, 190, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 500, 530,570, 580,610,  630, 660,710, 750, 780,800, 830, 850]
if options.width == "_w10.00":EXPmassesT = [ 150,160, 170, 180, 190, 200,220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 430, 460, 500, 530,560, 600,610,640, 660, 750, 770,800,  850]


OBSmassesT = [150,200, 250,300, 350,400, 450,500, 550,600,650,700, 750,800,  850]
EXPmassesT = [150,200, 250,300, 350,400, 450,500, 550,600,650,700, 750,800,  850]






#OBSmassesT = numpy.arange(110,150.1,1)
#EXPmassesT = numpy.arange(110,150.1,1)
epsilon = 0.001  # make this smaller than your smallest step size

for m in OBSmassesT:
        #if "%.1f"%m=="%d.0"%(m+epsilon):continue   # sigh!
    OBSmasses.append(m)
    EXPmasses.append(m)

# Plotting Styles --------------------------------------------------------
OFFSETLOW=0
OFFSETHIGH=0
FONTSIZE=0.035
FILLSTYLE=1001
SMFILLSTYLE=3244
FILLCOLOR_95=ROOT.kYellow
FILLCOLOR_68=ROOT.kGreen
RANGEYABS=[0.0001,0.5]
RANGEYRAT=[0.1,2000]
#RANGEYRAT=[0.0, 2.2]
RANGEMU = [-4,3.0]
MINPV = 0.5*10E-5
MAXPV = 1.0
Lines = [1.,2.,3.]
MINMH=int(min(EXPmasses))
MAXMH=int(max(EXPmasses))

if options.show : ROOT.gROOT.SetBatch(False)
if options.addline and not options.pval : sys.exit("Cannot addlines unless running in pvalue")

# ------------------------------------------------------------------------
# SM Signal Normalizer
#if not options.doRatio:
#    ROOT.gROOT.ProcessLine(".L ../libLoopAll.so")
#    signalNormalizer = ROOT.Normalization_8TeV()
extraString = "SM"
# ------------------------------------------------------------------------
if options.pval:
     EXPmasses=[]
     options.doRatio=True

if options.Method=="MaxLikelihoodFit": 
    options.doRatio = True

if not options.doRatio and options.Method != "Frequentist": 
    ROOT.gROOT.ProcessLine(".L medianCalc.C+g")
    from ROOT import medianCalc
    from ROOT import FrequentistLimits
    

if options.bayes:
  BayesianFile =ROOT.TFile(options.bayes) 
  bayesObs = BayesianFile.Get("observed")
 
Method = options.Method
if not options.path: options.path=Method


EXPName = options.path+"/higgsCombineEXPECTED."+Method
if Method =="MaxLikelihoodFit": EXPName = options.path+"/higgsCombineTest."+Method
if Method == "Asymptotic" or Method == "AsymptoticNew":  EXPName = options.path+"/higgsCombine"+options.width+"."+Method  # everyhting contained here
if Method == "ProfileLikelihood" or Method=="Asymptotic" or Method=="AsymptoticNew" or Method=="MaxLikelihoodFit":
  OBSName = options.path+"/higgsCombineTest."+Method
if Method == "Bayesian":
  OBSName = options.path+"/higgsCombineOBSERVED.MarkovChainMC"
if Method == "HybridNew":
  OBSName = options.path+"/higgsCombineOBSERVED.HybridNew"


if Method == "HybridNew" or Method == "Asymptotic" or Method == "AsymptoticNew": EXPmasses = OBSmasses[:]
# Done with setip
# --------------------------------------------------------------------------
ROOT.gROOT.ProcessLine( \
   "struct Entry{   \
    double r;       \
   };"
)
from ROOT import Entry
def getOBSERVED(file,entry=0):
  try:
   tree = file.Get("limit")
  except:
   return -1
  br = tree.GetBranch("limit")
  c = Entry()
  br.SetAddress(ROOT.AddressOf(c,'r'))
  tree.GetEntry(entry)
  return c.r

if Method=="HybridNew":
  EXPfiles=[]
  EXPmasses = OBSmasses[:]
  for m in EXPmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon):    # sigh!
      EXPfiles.append(ROOT.TFile(EXPName+".mH%d.quant0.500.root"%m))
    else:
      EXPfiles.append(ROOT.TFile(EXPName+".mH%.1f.quant0.500.root"%m))
    if options.verbose: print "expected MH - ", m, "File - ", EXPfiles[-1].GetName()
  
elif Method=="Asymptotic" or Method=="AsymptoticNew" or Method=="MaxLikelihoodFit":
  EXPfiles=[]
  EXPmasses = OBSmasses[:]
  for m in EXPmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon) and not options.sideband:   # sigh!
      EXPfiles.append(ROOT.TFile(EXPName+".mH%d.123456.root"%(m+epsilon)))
    else:
      EXPfiles.append(ROOT.TFile(EXPName+".mH%.1f.123456.root"%m))
    if options.verbose: print "expected MH - ", m, "File - ", EXPfiles[-1].GetName()

else:
  EXPfiles=[]
  for m in EXPmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon) and not options.sideband:   # sigh!
      EXPfiles.append(ROOT.TFile(EXPName+".mH%d.root"%(m+epsilon)))
    else:
      EXPfiles.append(ROOT.TFile(EXPName+".mH%.1f.root"%m))
    if options.verbose: print "expected MH - ", m, "File - ", EXPfiles[-1].GetName()

# Get the observed limits - Currently only does up to 1 decimal mass points
OBSfiles = []
if not options.expectedOnly:
  for m in OBSmasses:
    if "%.1f"%m=="%d.0"%(m+epsilon) and not options.sideband:   # sigh!
      OBSfiles.append(ROOT.TFile(OBSName+".mH%d.root"%(m+epsilon)))
    else:
      OBSfiles.append(ROOT.TFile(OBSName+".mH%.1f.root"%m))
    if options.verbose: print "observed MH - ", m, "File - ", OBSfiles[-1].GetName()

  if Method == "Asymptotic" or Method =="AsymptoticNew" :  obs = [getOBSERVED(O,5) for O in OBSfiles] # observed is last entry in these files
  else: obs = [getOBSERVED(O) for O in OBSfiles]

else:
  obs = [0 for O in OBSmasses]
  OBSfiles = obs[:]

# -------------------------------------------------------------------------------------------------------------------------------------------
# Set-up the GRAPHS

graph68  = ROOT.TGraphAsymmErrors()
graph95  = ROOT.TGraphAsymmErrors()
graphMed = ROOT.TGraphAsymmErrors()
graphObs = ROOT.TGraphAsymmErrors()
graphOne = ROOT.TGraphAsymmErrors()
dummyGraph= ROOT.TGraphAsymmErrors()

graph68up = ROOT.TGraphErrors()
graph68dn = ROOT.TGraphErrors()
graph95up = ROOT.TGraphErrors()
graph95dn = ROOT.TGraphErrors()
graphmede = ROOT.TGraphErrors()

graph68.SetLineColor(1)
graph95.SetLineColor(1)
graph68.SetLineStyle(2)
graph95.SetLineStyle(2)
graph68.SetLineWidth(2)
graph95.SetLineWidth(2)


MG = ROOT.TMultiGraph()

def MakeMlfPlot(MG):
    legend=ROOT.TLegend(0.55,0.65,0.89,0.85) #0.15-0.3
    legend.SetFillColor(10)
    legend.SetTextFont(42)
    legend.SetTextSize(FONTSIZE)
    graph68.SetLineStyle(1)
    legend.AddEntry(graph68,"#pm 1#sigma Uncertainty","F")

    if options.square : c = ROOT.TCanvas("c","c",600,600)
    else :c = ROOT.TCanvas("c","c",800,600)

    dhist = ROOT.TH1F("dh","dh",100,MINMH,MAXMH)
    dhist.GetYaxis().SetTitleOffset(1.2)
    dhist.GetXaxis().SetTitleOffset(1.2)
    dhist.GetYaxis().SetTitleSize(0.04)
    dhist.GetXaxis().SetTitleSize(0.04)
    dhist.GetYaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetRangeUser(MINMH,MAXMH)
    dhist.GetYaxis().SetRangeUser(RANGEMU[0],RANGEMU[1])
    dhist.GetXaxis().SetTitle("m_{H} (GeV)")
    dhist.GetYaxis().SetTitle("Best fit #sigma/#sigma_{SM}")
    dhist.Draw("AXIS")

    MG.Draw("L3")

    # ------------------------------------------------------------------------
    # Additional Lines stored in --addline -----------------------------------
    for lineF in options.addline:

        # Parse the string, should be file.root:color:linestyle:legend entry    
        vals = lineF.split(":")
        ftmp = ROOT.TFile(vals[0])
        grext = ftmp.Get("observed")
        grext.SetLineColor(int(vals[1]))
        grext.SetLineStyle(int(vals[2]))
        grext.SetLineWidth(2)
        legend.AddEntry(grext,vals[3],"L")
        grext.Draw("same")
    # ------------------------------------------------------------------------
        

    c.Update()
    text = ROOT.TLatex()
    text.SetTextColor(ROOT.kRed)
    text.SetTextSize(FONTSIZE)
    text.SetTextFont(42)

    if options.doRatio: graphOne.Draw("L")
    c.SetGrid(not options.nogrid)
    if not options.nogrid: dhist.Draw("AXIGSAME")

    mytext= ROOT.TLatex()
    mytext.SetTextSize(FONTSIZE)
    mytext.SetTextFont(42)
    mytext.SetNDC()

    mytext.DrawLatex(0.18,0.24,"CMS Preliminary")
    for t,lineT in enumerate(options.addtxt):
        mytext.DrawLatex(0.18,0.23-(t+1)*0.04,"%s"%(lineT))
    legend.Draw()
    ROOT.gPad.RedrawAxis();
    
    if options.show:raw_input("Looks Ok?")
    c.SaveAs("maxlhplot.pdf")
    c.SaveAs("maxlhplot.png")

#-------------------------------------------------------------------------
def MakePvalPlot(MG):

    legend=ROOT.TLegend(0.55,0.65,0.89,0.85) #0.15-0.3
    legend.SetFillColor(10)
    legend.SetTextFont(42)
    legend.SetTextSize(FONTSIZE)
    legend.AddEntry(graphObs,"Observed","L")

    if options.square : c = ROOT.TCanvas("c","c",600,600)
    else :c = ROOT.TCanvas("c","c",800,600)

    dhist = ROOT.TH1F("dh","dh",100,MINMH,MAXMH)
    dhist.GetYaxis().SetTitleOffset(1.5)
    dhist.GetXaxis().SetTitleOffset(1.2)
    dhist.GetYaxis().SetTitleSize(0.04)
    dhist.GetXaxis().SetTitleSize(0.04)
    dhist.GetYaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetRangeUser(MINMH,MAXMH)
    dhist.GetYaxis().SetRangeUser(MINPV,MAXPV)
    dhist.GetXaxis().SetTitle("m_{H} (GeV)")
    dhist.GetYaxis().SetTitle("Local p-value")
    dhist.Draw("AXIS")

    MG.Draw("L3")

    # ------------------------------------------------------------------------
    # Additional Lines stored in --addline -----------------------------------
    for lineF in options.addline:

        # Parse the string, should be file.root:color:linestyle:legend entry    
        vals = lineF.split(":")
        ftmp = ROOT.TFile(vals[0])
        grext = ftmp.Get("observed")
        grext.SetLineColor(int(vals[1]))
        grext.SetLineStyle(int(vals[2]))
        grext.SetLineWidth(2)
        legend.AddEntry(grext,vals[3],"L")
        grext.Draw("same")
    # ------------------------------------------------------------------------
        

    c.Update()
    text = ROOT.TLatex()
    text.SetTextColor(ROOT.kRed)
    text.SetTextSize(FONTSIZE)
    text.SetTextFont(42)
    
        
    Vals=[ROOT.RooStats.SignificanceToPValue(L) for L in Lines]
    TLines=[ROOT.TLine(MINMH,V,MAXMH,V) for V in Vals]

    for j,TL in enumerate(TLines):
        TL.SetLineStyle(1)
        TL.SetLineColor(2)
        TL.SetLineWidth(1)
        TL.Draw("same")
        text.DrawLatex(MAXMH+0.2,Vals[j]*0.88,"%d #sigma"%Lines[j])

    c.SetGrid(not options.nogrid)
    c.SetLogy();
    if not options.nogrid: dhist.Draw("AXIGSAME")

    mytext= ROOT.TLatex()
    mytext.SetTextSize(FONTSIZE)
    mytext.SetTextFont(42)
    mytext.SetNDC()

    mytext.DrawLatex(0.18,0.24,"CMS Preliminary")
    for t,lineT in enumerate(options.addtxt):
        mytext.DrawLatex(0.18,0.23-(t+1)*0.04,"%s"%(lineT))
    legend.Draw()
    ROOT.gPad.RedrawAxis();
    
    if options.show:raw_input("Looks Ok?")
    c.SaveAs("pvaluesplot.pdf")
    c.SaveAs("pvaluesplot.png")
#-------------------------------------------------------------------------

def MakeLimitPlot(MG):
    cat_str="Cat: "
    leg=ROOT.TLegend(0.5,0.65,0.89,0.85) #0.2 0.2-0.59 0.4
    if options.cat=="CH0": cat_str+="0"
    if options.cat=="CH1": cat_str+="1"
    if options.cat=="CH2": cat_str+="2"
    if options.cat=="CH3": cat_str+="3"
    if options.cat=="COMB": cat_str="All Events"
    if options.cat=="COMBWSyst": cat_str="All Events"
    if options.cat=="COMBOnlyEB": cat_str="Events in the Barrel"
    
    if options.width == "_w0.10":leg.SetHeader("#splitline{"+cat_str+"}{Width: 0.1 [GeV]}")
    if options.width == "_w2.00":leg.SetHeader("#splitline{"+cat_str+"}{Width: 2% of m_{H} [GeV]}")
    if options.width == "_w5.00":leg.SetHeader("#splitline{"+cat_str+"}{Width: 5% of m_{H} [GeV]}")
    if options.width == "_w7.00":leg.SetHeader("#splitline{"+cat_str+"}{Width: 7% of m_{H} [GeV]}")
    if options.width == "_w10.00":leg.SetHeader("#splitline{"+cat_str+"}{Width: 10% of m_{H} [GeV]}")
    
    leg.SetFillColor(10)

    # Different entries for the different methods
    LegendEntry = ""
    if Method == "ProfileLikelihood": LegendEntry = "PL"
    if Method == "Bayesian": LegendEntry = "Bayesian"
    if Method == "HybridNew": LegendEntry = "CLs"
    if Method == "Asymptotic": LegendEntry = "CLs (Asymptotic)"
    if Method == "AsymptoticNew": LegendEntry = "CLs (Asymptotic)"

    if not options.expectedOnly: leg.AddEntry(graphObs,"Observed","L")
    if options.bayes and not options.expectedOnly: leg.AddEntry(bayesObs,"Observed Bayesian Limit","L")
    leg.AddEntry(graph68,"Expected #pm 1#sigma","FL")
    leg.AddEntry(graph95,"Expected #pm 2#sigma","FL")

    leg.SetTextFont(42)
    leg.SetTextSize(FONTSIZE)

    if options.square : C = ROOT.TCanvas("c","c",600,600)
    else: C = ROOT.TCanvas("c","c",700,600)

    C.SetGrid(not options.nogrid)
    dummyHist = ROOT.TH1D("dummy","",1,min(OBSmasses)-OFFSETLOW,max(OBSmasses)+OFFSETHIGH)
    dummyHist.SetTitleSize(0.04,"XY")
    dummyHist.Draw("AXIS")
    MG.Draw("L3")
    if not options.nogrid: dummyHist.Draw("AXIGSAME")

    dummyHist.GetXaxis().SetTitle("m_{H} [GeV]")
    dummyHist.GetXaxis().SetRangeUser(min(OBSmasses)-OFFSETLOW,max(OBSmasses)+OFFSETHIGH)
    if options.doRatio:
     dummyHist.GetYaxis().SetRangeUser(RANGEYRAT[0],RANGEYRAT[1])
     dummyHist.GetYaxis().SetNdivisions(5,int("%d"%(RANGEYRAT[1]-RANGEYRAT[0])),0)
     dummyHist.GetYaxis().SetTitle("\sigma(H#rightarrow #gamma #gamma)_{95%%CL} / \sigma(H#rightarrow #gamma #gamma)_{%s}"%extraString)
    else: 
     dummyHist.GetYaxis().SetRangeUser(RANGEYABS[0],RANGEYABS[1])
     dummyHist.GetYaxis().SetNdivisions(5,int("%d"%(RANGEYABS[1]-RANGEYABS[0])),0)
     dummyHist.GetYaxis().SetTitle("\sigma #times BR(H#rightarrow #gamma #gamma)_{95%CL} (pb)")

    dummyHist.GetYaxis().SetTitleOffset(1.3)
    dummyHist.GetXaxis().SetTitleOffset(1.25)

    MG.SetTitle("")
    mytext = ROOT.TLatex()
    mytext.SetTextSize(FONTSIZE)

    mytext.SetNDC()
    mytext.SetTextFont(42)
    mytext.SetTextSize(FONTSIZE)
    mytext.DrawLatex(0.12,0.96,"CMS Preliminary, 19.5 fb^{-1}                                                   #sqrt{s} = 8 TeV")
    for t,lineT in enumerate(options.addtxt):
        mytext.DrawLatex(0.16,0.84-(t+1)*(0.04),"%s"%lineT)
  
    leg.Draw()
    ROOT.gPad.RedrawAxis();
    if options.show:raw_input("Looks Ok?")
    C.SetLogy()
    #Make a bunch of extensions to the plots
    outputname=options.path+"/limit"
    www_outputname=www_path+"/limit"
    if options.doSmooth: outputname+="_smooth"
    if options.doSmooth: www_outputname+="_smooth"
    outputname+="_"+Method
    www_outputname+="_"+Method
    if options.doRatio: outputname+="_ratio"
    if options.doRatio: www_outputname+="_ratio"
    if options.append!="": outputname+="_"+options.append
    if options.append!="": www_outputname+="_"+options.append
    if options.width == "_w0.10": outputname+="_w01"
    if options.width == "_w0.10": www_outputname+="_w01"
    if options.width == "_w2.00": outputname+="_w2"
    if options.width == "_w2.00": www_outputname+="_w2"
    if options.width == "_w5.00": outputname+="_w5"
    if options.width == "_w5.00": www_outputname+="_w5"
    if options.width == "_w7.00": outputname+="_w7"
    if options.width == "_w7.00": www_outputname+="_w7"
    if options.width == "_w10.00": outputname+="_w10"
    if options.width == "_w10.00": www_outputname+="_w10"
    
    
    outputname+="_"+options.cat
    www_outputname+="_"+options.cat
    types=[".pdf",".png",".eps"]
    for type in types: C.SaveAs(outputname+type)
    for type in types: C.SaveAs(www_outputname+type)

#-------------------------------------------------------------------------


#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#EXPECTED + Bands
for i,mass,f in zip(range(len(EXPfiles)),EXPmasses,EXPfiles):
  if options.pval: continue
  sm = 1.
  xsec150=(13.65+1.280+0.3681+0.2159+0.07403)
  median = array.array('d',[0])
  up68   = array.array('d',[0])
  dn68   = array.array('d',[0])
  up95   = array.array('d',[0])
  dn95   = array.array('d',[0])

  #if not options.doRatio:
  if mass == 110:
      sm=0.00195*(25.04+1.809+1.060+0.5869+0.1887)
      sm150=0.00195*xsec150
  if mass == 120:
      sm=0.00223*(21.13+1.649+0.7966+0.4483+0.1470)
      sm150=0.00223*xsec150
  if mass == 130:
      sm=0.00225*(18.07+1.511+0.6095+0.3473+0.1157)
      sm150=0.00225*xsec150
  if mass==140:
      sm=0.00193*(15.63+1.389+0.4713+0.2728+0.09207)
      sm150=0.00193*xsec150
  if mass == 150:
      sm = 0.00137*(13.65+1.280+0.3681+0.2159+0.07403)
      sm150=0.00137*xsec150
  if mass == 155:
      sm = 0.000998*(12.79+1.231+0.3252+0.1923+0.06664)
      sm150 = 0.000998*xsec150
  if mass == 160:
      sm = 0.000532*(11.95+1.185+0.2817+0.1687+0.06013)
      sm150 = 0.000532*xsec150
  if mass == 165:
      sm = 0.000229*(10.89+1.141+0.2592+0.1561+0.05439)
      sm150 = 0.000229*xsec150
  if mass == 170:
      sm = 0.000158*(10.12+1.098+0.2329+0.1408+0.04930)
      sm150 = 0.000158*xsec150
  if mass == 175:
      sm = 0.000127*(9.476+1.055+0.2089+0.1266+0.04480)
      sm150 = 0.000127*xsec150
  if mass == 180:
      sm = 0.000105*(8.874+1.015+0.1883+0.1137+0.04080)
      sm150 = 0.000105*xsec150   
  if mass == 185:
      sm = 0.0000843*(8.326+0.9760+0.1715+0.1038+0.03725)
      sm150= 0.0000843*xsec150
  if mass == 195:
      sm = 0.0000617*(7.868+0.9018+0.1416+0.08584+0.03125)
      sm150 = 0.0000617*xsec150
  if mass == 200:
      sm = 5.5e-05*(7.127+0.8685+0.1286+0.07827+0.02872)
      sm150 = 5.5e-05*xsec150
  if mass == 210:
      sm = 4.55e-05*(6.534+0.8163+0.1070+0.06526+0.02442)
      sm150 = 4.55e-05*xsec150
  if mass == 220:
      sm = 3.83e-05*(6.038+0.7677+0.08961+0.05476+0.02094)
      sm150 = 3.83e-05*xsec150
  if mass == 230:
      sm = 3.28e-05*(5.593+0.7190+0.07559+0.04614+0.01810)
      sm150 = 3.28e-05*xsec150
  if mass == 240:
      sm = 2.84e-05*(5.183+0.6703+0.06398+0.03901+0.01574)
      sm150 = 2.84e-05*xsec150
  if mass == 250:
      sm = 2.47e-05*(4.802+0.6225+0.05450+0.03314+0.01380)
      sm150 = 2.47e-05*xsec150
  if mass == 260:
      sm = 2.17e-05*(4.479+0.5797+0.04660+0.02821+0.01219)
      sm150 = 2.17e-05*xsec150
  if mass == 270:
      sm = 1.91e-05*(4.198+0.5401+0.04016+0.02415+0.01083)
      sm150 = 1.91e-05*xsec150
  if mass == 280:
      sm = 1.68e-05*(3.964+0.5045+0.03459+0.02075+0.009686)
      sm150 = 1.68e-05*xsec150
  if mass == 290:
      sm = 1.5e-05*(3.767+0.4716+0.02993+0.01788+0.008705)
      sm150 = 1.5e-05*xsec150
  if mass == 300:
      sm = 1.3e-05*(3.606+0.4408+0.02602+0.01547+0.007862)
      sm150 = 1.3e-05*xsec150
  if mass == 310:
      sm = 1.2e-05*(3.482+0.4131)
      sm150 = 1.2e-05*xsec150
  if mass == 320:
      sm = 1.09e-05*(3.392+0.3875)
      sm150 = 1.09e-05*xsec150
  if mass == 330:
      sm = 9.83e-06*(3.349+0.3637)
      sm150 = 9.83e-06*xsec150
  if mass == 340:
      sm = 8.92e-06*(3.367+0.3422)
      sm150 = 8.92e-06*xsec150
  if mass == 350:
      sm = 7.73e-06*(3.406+0.3200)
      sm150 = 7.73e-06*xsec150
  if mass == 360:
      sm = 6.25e-06*(3.390+0.3028)
      sm150 = 6.25e-06*xsec150
  if mass == 370:
      sm = 5.04e-06*(3.336+0.2896)
      sm150 = 5.04e-06*xsec150
  if mass == 380:
      sm = 4.06e-06*(3.235+0.2776)
      sm150 = 4.06e-06*xsec150
  if mass == 390:
      sm = 3.28e-06*(3.093+0.2660)
      sm150 = 3.28e-06*xsec150
  if mass == 400:
      sm = 2.65e-06*(2.924+0.2543)
      sm150 = 2.65e-06*xsec150
  if mass == 430:
      sm = 0.00000141*(2.3+0.22)
      sm150 = 0.00000141*xsec150
  if mass == 460:
      sm = 0.000000751*(1.837+0.1905)
      sm150 =0.000000751*xsec150
  if mass == 500:
      sm = 0.000000312*(1.28+0.1561)
      sm150 =0.000000312*xsec150
  if mass == 530:
      sm = 0.000000156*(0.85+0.13)
      sm150 =0.000000156*xsec150
  if mass == 550:
      sm = 0.0000000990*(0.8141+0.1223)
      sm150 = 0.0000000990*xsec150
  if mass == 560:
      sm = 0.0000000810*(0.744+0.1166)
      sm150 =0.0000000810*xsec150
  if mass == 600:
      sm = 0.0000000503*(0.523+0.09688)
      sm150 =0.0000000503*xsec150
  if mass == 630:
      sm = 0.0000000493*(0.4+0.085)
      sm150 =0.0000000493*xsec150
  if mass == 650:
      sm = 0.0000000582*(0.3423+0.07784)
      sm150 = 0.0000000582*xsec150
  if mass == 660:
      sm = 0.0000000660*(0.3152+0.07459)
      sm150 = 0.0000000660*xsec150
  if mass == 700:
      sm = 0.000000111*(0.2288+0.06330)
      sm150 = 0.000000111*xsec150
  if mass == 730:
      sm = 0.000000152*(0.18+0.056)
      sm150 =0.000000152*xsec150
  if mass == 760:
      sm = 0.000000192*(0.1455+0.05032)
      sm150 =0.000000192*xsec150
  if mass == 800:
      sm = 0.000000244*(0.1095+0.04365)
      sm150 =0.000000244*xsec150
  if mass == 830:
      sm = 0.000000282*(0.09+0.039)
      sm150 =0.000000282*xsec150
  if mass == 850:
      sm = 0.000000315*(0.07811+0.037)
      sm150 =0.000000315*xsec150
  sm150=(13.65+1.280+0.3681+0.2159+0.07403)*0.00137
  if options.doRatio: sm150=sm150/sm
  
  #if not options.doRatio: sm = signalNormalizer.GetBR(mass)*signalNormalizer.GetXsection(mass)
  if Method == "Asymptotic" or Method=="AsymptoticNew":   
      median[0] = getOBSERVED(f,2)
      up95[0]   = getOBSERVED(f,4)
      dn95[0]   = getOBSERVED(f,0)
      up68[0]   = getOBSERVED(f,3)
      dn68[0]   = getOBSERVED(f,1)

  elif Method=="MaxLikelihoodFit":
      median[0] = getOBSERVED(f,0)
      up95[0]   = median[0]
      dn95[0]   = median[0]
      up68[0]   = getOBSERVED(f,2)
      dn68[0]   = getOBSERVED(f,1)

  else:
    tree = f.Get("limit")
    medianCalc("r_mH"+str(mass),tree,median,up68,dn68,up95,dn95)

  graph68.SetPoint(i,float(mass),median[0]*sm150)
  graph95.SetPoint(i,float(mass),median[0]*sm150)
  graphMed.SetPoint(i,float(mass),median[0]*sm150)
  
  print str(mass)+"    "+str(median[0]*sm150)+"   "+str(median[0]*sm150*19620.)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+" "+options.width
  
  if options.width == "_w0.1" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   0.1\n")
  if options.width == "_w0.10" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   0.1\n")
  if options.width == "_w2.00" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   2\n")
  if options.width == "_w5.00" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   5\n")
  if options.width == "_w7.00" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   7\n")
  if options.width == "_w10.00" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   10\n")
  
  if options.doRatio: graphOne.SetPoint(i,float(mass),1.*sm150)
  
  if Method == "HybridNew":

      up95[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.975.root"))
      dn95[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.025.root"))
      up68[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.840.root"))
      dn68[0]   = FrequentistLimits(f.GetName().replace("0.500.root","0.160.root"))

  
  diff95_up = abs(median[0] - up95[0])*sm150
  diff95_dn = abs(median[0] - dn95[0])*sm150
  diff68_up = abs(median[0] - up68[0])*sm150
  diff68_dn = abs(median[0] - dn68[0])*sm150
  
  graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
  graph95.SetPointError(i,0,0,diff95_dn,diff95_up)
  graphMed.SetPointError(i,0,0,0,0)
  graphOne.SetPointError(i,0,0,0,0)

  if options.doSmooth:  # Always fit the absolute not the ratio
   # sm=1.
    print str(mass)+"    "+str(median[0]*sm150)+"   "+options.width
    graphmede.SetPoint(i,float(mass),median[0]*sm150)
    graph68up.SetPoint(i,float(mass),up68[0]*sm150)
    graph68dn.SetPoint(i,float(mass),dn68[0]*sm150)
    graph95up.SetPoint(i,float(mass),up95[0]*sm150)
    graph95dn.SetPoint(i,float(mass),dn95[0]*sm150)

# Smooth the Bands set -doSmooth
# Since i always fitted to the Absolute, need to see if i want the Ratio instead
if options.doSmooth:
  print str(mass)+"    "+str(median[0]*sm150)+"   "+options.width
  fitstring = "[0] + [1]*x*x + [2]*x*x*x +[3]*x*x*x*x + [4]*x"
  medfunc = ROOT.TF1("medfunc",fitstring,  150.,850.);
  up68func = ROOT.TF1("up68func",fitstring,150.,850.);
  dn68func = ROOT.TF1("dn68func",fitstring,150.,850.);
  up95func = ROOT.TF1("up95func",fitstring,150.,850.);
  dn95func = ROOT.TF1("dn95func",fitstring,150.,850.);

  graphmede.Fit(medfunc,"R,M,EX0","Q")
  graph68up.Fit(up68func,"R,M,EX0","Q")
  graph68dn.Fit(dn68func,"R,M,EX0","Q")
  graph95up.Fit(up95func,"R,M,EX0","Q")
  graph95dn.Fit(dn95func,"R,M,EX0","Q")
  
  newCanvas = ROOT.TCanvas()
  graphmede.SetMarkerSize(0.8)
  graphmede.SetMarkerStyle(21)
  graphmede.Draw("A")
  newCanvas.SaveAs("smoothTest.pdf")
    
  for i,mass in zip(range(len(EXPmasses)),EXPmasses):
     print str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+" "+options.width 
     sm=1.0
     sm150=1.0
  
  #if not options.doRatio:
  # sm = signalNormalizer.GetBR(mass)*signalNormalizer.GetXsection(mass)
 
  #if not options.doRatio:
     if mass == 110:
          sm=0.00195*(25.04+1.809+1.060+0.5869+0.1887)
          sm150=0.00195*xsec150
     if mass == 120:
          sm=0.00223*(21.13+1.649+0.7966+0.4483+0.1470)
          sm150=0.00223*xsec150
     if mass == 130:
          sm=0.00225*(18.07+1.511+0.6095+0.3473+0.1157)
          sm150=0.00225*xsec150
     if mass==140:
          sm=0.00193*(15.63+1.389+0.4713+0.2728+0.09207)
          sm150=0.00193*xsec150
     if mass == 150:
          sm = 0.00137*(13.65+1.280+0.3681+0.2159+0.07403)
          sm150=0.00137*xsec150
     if mass == 155:
          sm = 0.000998*(12.79+1.231+0.3252+0.1923+0.06664)
          sm150 = 0.000998*xsec150
     if mass == 160:
          sm = 0.000532*(11.95+1.185+0.2817+0.1687+0.06013)
          sm150 = 0.000532*xsec150
     if mass == 165:
          sm = 0.000229*(10.89+1.141+0.2592+0.1561+0.05439)
          sm150 = 0.000229*xsec150
     if mass == 170:
          sm = 0.000158*(10.12+1.098+0.2329+0.1408+0.04930)
          sm150 = 0.000158*xsec150
     if mass == 175:
          sm = 0.000127*(9.476+1.055+0.2089+0.1266+0.04480)
          sm150 = 0.000127*xsec150
     if mass == 180:
          sm = 0.000105*(8.874+1.015+0.1883+0.1137+0.04080)
          sm150 = 0.000105*xsec150   
     if mass == 185:
          sm = 0.0000843*(8.326+0.9760+0.1715+0.1038+0.03725)
          sm150= 0.0000843*xsec150
     if mass == 195:
          sm = 0.0000617*(7.868+0.9018+0.1416+0.08584+0.03125)
          sm150 = 0.0000617*xsec150
     if mass == 200:
          sm = 5.5e-05*(7.127+0.8685+0.1286+0.07827+0.02872)
          sm150 = 5.5e-05*xsec150
     if mass == 210:
          sm = 4.55e-05*(6.534+0.8163+0.1070+0.06526+0.02442)
          sm150 = 4.55e-05*xsec150
     if mass == 220:
          sm = 3.83e-05*(6.038+0.7677+0.08961+0.05476+0.02094)
          sm150 = 3.83e-05*xsec150
     if mass == 230:
          sm = 3.28e-05*(5.593+0.7190+0.07559+0.04614+0.01810)
          sm150 = 3.28e-05*xsec150
     if mass == 240:
          sm = 2.84e-05*(5.183+0.6703+0.06398+0.03901+0.01574)
          sm150 = 2.84e-05*xsec150
     if mass == 250:
          sm = 2.47e-05*(4.802+0.6225+0.05450+0.03314+0.01380)
          sm150 = 2.47e-05*xsec150
     if mass == 260:
          sm = 2.17e-05*(4.479+0.5797+0.04660+0.02821+0.01219)
          sm150 = 2.17e-05*xsec150
     if mass == 270:
          sm = 1.91e-05*(4.198+0.5401+0.04016+0.02415+0.01083)
          sm150 = 1.91e-05*xsec150
     if mass == 280:
          sm = 1.68e-05*(3.964+0.5045+0.03459+0.02075+0.009686)
          sm150 = 1.68e-05*xsec150
     if mass == 290:
          sm = 1.5e-05*(3.767+0.4716+0.02993+0.01788+0.008705)
          sm150 = 1.5e-05*xsec150
     if mass == 300:
          sm = 1.3e-05*(3.606+0.4408+0.02602+0.01547+0.007862)
          sm150 = 1.3e-05*xsec150
     if mass == 310:
          sm = 1.2e-05*(3.482+0.4131)
          sm150 = 1.2e-05*xsec150
     if mass == 320:
          sm = 1.09e-05*(3.392+0.3875)
          sm150 = 1.09e-05*xsec150
     if mass == 330:
          sm = 9.83e-06*(3.349+0.3637)
          sm150 = 9.83e-06*xsec150
     if mass == 340:
          sm = 8.92e-06*(3.367+0.3422)
          sm150 = 8.92e-06*xsec150
     if mass == 350:
          sm = 7.73e-06*(3.406+0.3200)
          sm150 = 7.73e-06*xsec150
     if mass == 360:
          sm = 6.25e-06*(3.390+0.3028)
          sm150 = 6.25e-06*xsec150
     if mass == 370:
          sm = 5.04e-06*(3.336+0.2896)
          sm150 = 5.04e-06*xsec150
     if mass == 380:
          sm = 4.06e-06*(3.235+0.2776)
          sm150 = 4.06e-06*xsec150
     if mass == 390:
          sm = 3.28e-06*(3.093+0.2660)
          sm150 = 3.28e-06*xsec150
     if mass == 400:
          sm = 2.65e-06*(2.924+0.2543)
          sm150 = 2.65e-06*xsec150
     if mass == 430:
         sm = 0.00000141*(2.3+0.22)
         sm150 = 0.00000141*xsec150
     if mass == 460:
         sm = 0.000000751*(1.837+0.1905)
         sm150 =0.000000751*xsec150
     if mass == 500:
         sm = 0.000000312*(1.28+0.1561)
         sm150 =0.000000312*xsec150
     if mass == 530:
         sm = 0.000000156*(0.85+0.13)
         sm150 =0.000000156*xsec150
     if mass == 550:
         sm = 0.0000000990*(0.8141+0.1223)
         sm150 = 0.0000000990*xsec150
     if mass == 560:
         sm = 0.0000000810*(0.744+0.1166)
         sm150 =0.0000000810*xsec150
     if mass == 600:
         sm = 0.0000000503*(0.523+0.09688)
         sm150 =0.0000000503*xsec150
     if mass == 630:
         sm = 0.0000000493*(0.4+0.085)
         sm150 =0.0000000493*xsec150
     if mass == 650:
         sm = 0.0000000582*(0.3423+0.07784)
         sm150 = 0.0000000582*xsec150
     if mass == 660:
         sm = 0.0000000660*(0.3152+0.07459)
         sm150 = 0.0000000660*xsec150
     if mass == 700:
         sm = 0.000000111*(0.2288+0.06330)
         sm150 = 0.000000111*xsec150
     if mass == 730:
         sm = 0.000000152*(0.18+0.056)
         sm150 =0.000000152*xsec150
     if mass == 760:
         sm = 0.000000192*(0.1455+0.05032)
         sm150 =0.000000192*xsec150
     if mass == 800:
         sm = 0.000000244*(0.1095+0.04365)
         sm150 =0.000000244*xsec150
     if mass == 830:
         sm = 0.000000282*(0.09+0.039)
         sm150 =0.000000282*xsec150
     if mass == 850:
         sm = 0.000000315*(0.07811+0.037)
         sm150 =0.000000315*xsec150
     sm150=(13.65+1.280+0.3681+0.2159+0.07403)*0.00137    
     if options.doRatio: sm150=sm150/sm
     print str(mass)+"    "+str(median[0]*sm150)+"   "+str(median[0]*sm150*19620.)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+" "+options.width
     
     if options.width == "_w0.1" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   0.1\n")
     if options.width == "_w0.10" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   0.1\n")
     if options.width == "_w2.00" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   2\n")
     if options.width == "_w5.00" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   5\n")
     if options.width == "_w7.00" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   7\n")
     if options.width == "_w10.00" : out_file.write(str(mass)+"    "+str(median[0]*sm150)+"  "+str(up68[0]*sm150)+"  "+str(dn68[0]*sm150)+"   10\n")
 
     mediansmooth = medfunc.Eval(mass)
      
     graphMed.SetPoint(i,mass,mediansmooth*sm150)
     graph68.SetPoint(i,mass,mediansmooth*sm150)
     graph95.SetPoint(i,mass,mediansmooth*sm150)
     print mediansmooth*sm150
     diff95_up = abs(mediansmooth - up95func.Eval(mass))*sm150
     diff95_dn = abs(mediansmooth - dn95func.Eval(mass))*sm150
     diff68_up = abs(mediansmooth - up68func.Eval(mass))*sm150
     diff68_dn = abs(mediansmooth - dn68func.Eval(mass))*sm150
      
     graph68.SetPointError(i,0,0,diff68_dn,diff68_up)
     graph95.SetPointError(i,0,0,diff95_dn,diff95_up)

# OBSERVED -------- easy as that !
for i,mass in zip(range(len(OBSfiles)),OBSmasses):

    sm = 1.;
    if obs[i] ==-1: continue
    if not options.doRatio: sm = 1.#signalNormalizer.GetBR(M)*signalNormalizer.GetXsection(M)
    graphObs.SetPoint(i,float(mass),obs[i]*sm150)
    graphObs.SetPointError(i,0,0,0,0)

# Finally setup the graphs and plot
graph95.SetFillColor(FILLCOLOR_95)
graph95.SetFillStyle(FILLSTYLE)
graph68.SetFillColor(FILLCOLOR_68)
graph68.SetFillStyle(FILLSTYLE)
graphMed.SetLineStyle(2)
graphMed.SetLineColor(2)
graphMed.SetMarkerColor(2)
graphMed.SetLineWidth(3)
graphObs.SetLineWidth(3)

if options.bayes:
 bayesObs.SetLineWidth(3)
 bayesObs.SetLineColor(4)
 bayesObs.SetMarkerColor(4)
 bayesObs.SetLineStyle(7)

graphOne.SetLineWidth(3)
graphOne.SetLineColor(ROOT.kRed)
graphOne.SetMarkerColor(ROOT.kRed)
graphObs.SetMarkerStyle(20)
graphObs.SetMarkerSize(2.0)
graphObs.SetLineColor(1)

graphMed.SetLineStyle(2)
graphMed.SetLineColor(ROOT.kBlack)
if not options.pval:MG.Add(graph95)
if not options.pval:MG.Add(graph68)
if not options.pval:MG.Add(graphMed)

if not options.expectedOnly:
  MG.Add(graphObs)
  if options.bayes:
   MG.Add(bayesObs)

#if not options.pval: MG.Add(graphOne) #livia

# Plot -------------------------------------
if options.pval: MakePvalPlot(MG)
elif Method=="MaxLikelihoodFit":  MakeMlfPlot(MG)
else:MakeLimitPlot(MG)
# ------------------------------------------
if options.outputLimits:
  print "Writing Limits To ROOT file --> ",options.outputLimits
  OUTTgraphs = ROOT.TFile(options.outputLimits,"RECREATE")
  graphObs.SetName("observed")
  graphObs.Write()

  if not options.pval:
   graphMed.SetName("median")
   graphMed.Write()
   graph68.SetName("sig1")
   graph68.Write()
   graph95.SetName("sig2")
   graph95.Write()

  OUTTgraphs.Write()
  out_file.close()
