#!/bin/bash

for mass in 450 #200 250 300 400 
  do
  for width in 10 #5 0
    do
    for min in 130 200 250 300 
      do
      for max in 600 700 800 1000
	do
   cat>>mk_workspace.C<<EOF
{
gROOT->ProcessLine(".L ProduceWorkspaces_Bias.C");
ProduceWorkspaces_Bias($mass, $width, $min, $max);
}
EOF

echo $mass $width $min $max
root -l -b -q mk_workspace.C 
     
      done # max
    done #min
  done #width
done #mass


