#!/bin/bash

for mass in 150 250 350 450 550 650 750 850 
  do
  for width in 0 10 #5 0
    do
    cat>>mk_workspace.C<<EOF

{
gROOT->ProcessLine(".L ProduceWorkspaces_Bias.C");
ProduceWorkspaces_Bias($mass, $width, 130, 1000);
}


EOF
    
    echo $mass $width 
    root -l -b -q mk_workspace.C > log_ws_$mass_$width.log
    rm mk_workspace.C 

  done #width
done #mass


