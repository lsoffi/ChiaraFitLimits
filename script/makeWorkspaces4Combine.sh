#!/bin/bash

for mass in 800 #200 300 400 500 600 700  
  do
  for width in 0.1 10 #5 0
    do
    cat>>mk_workspace.C<<EOF

{
gROOT->ProcessLine(".L ProduceWorkspaces.C");
ProduceWorkspaces($mass, $width);
}


EOF
    
    echo $mass $width 
    root -l -b -q mk_workspace.C > log_ws_$mass_$width.log
    rm mk_workspace.C 

  done #width
done #mass


