#!/bin/bash
#
#$ -N demHand1
#$ -m be
#$ -l qp=LOW,h_rt=24:00:00
#$ -cwd
nohup /usr/local/bin/matlab << EOF
HOME = getenv('HOME');
addpath([HOME '/mlprojects/matlab/general']);
cd ~/mlprojects/fgplvm/matlab
fgplvmToolboxes
demHand1
EOF
