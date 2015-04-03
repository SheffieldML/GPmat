#!/bin/bash
#
#$ -N demCmu35gplvm3
#$ -m be
#$ -l qp=LOW,h_rt=24:00:00
#$ -cwd
/usr/local/bin/matlab << EOF
HOME = getenv('HOME');
addpath([HOME '/mlprojects/matlab/general']);
cd ~/mlprojects/fgplvm/matlab
fgplvmToolboxes
demCmu35gplvm3
EOF
