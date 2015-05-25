#!/bin/bash

# Script to launch many PPA jobs at a PBS managed cluster.

function createMatlabFile
{
  # create matlab file
    cat > $1 <<EOF 
importTool('optimi', 0.12)
importTool('kern')
importTool('noise', 0.12)
importTool('ndlutil') 
importTool('ivm')
importTool('prior', 0.12)
%cd ~/mlprojects/ppa/matlab 
EOF
echo "ivmGunnarData($2, $3, $4, $5, $6);" >> $1
echo exit >> $1
}	 


function createScriptFile
{
  # create matlab file
    cat > $1 <<EOF 
#!/bin/bash
#
#PBS -m a
#PBS -M neil@ivy
#PBS -k oe
#PBS -l nodes=1:ppn=1
#PBS -q fast

# Send mail to me if job aborts -m a, -M neil@ivy
# Kill job if it is longer than an hour -l walltime=
EOF
echo "#PBS -N $2" >> $1
echo "" >> $1
echo "echo Running on " >> $1
echo "date" >> $1
echo "echo Using commands  " >> $1
echo "cat $3" >> $1
echo "unset DISPLAY" >> $1
echo "/usr/local/bin/matlab < $3" >> $1
}	 

NAME=ivmGunnarDataRBF
RUNNAME=$NAME
KERN="{'rbf', 'white', 'bias'}"

DATASETS=thyroid,100:titanic,100:heart,100:breast-cancer,100:banana,100:ringnorm,100:twonorm,100:waveform,100:diabetis,200:flare-solar,400:german,400:spice,400:image,400
DATANUM=`seq -s: 1 5`
INVWIDTHPARAMS=0.1:1:10

NOISE=\'probit\'
IFS=:
echo Creating run scripts ...
INDEX=$((0))
for dataNum in $DATANUM
  do
  for datasetDval in $DATASETS
    do
    for invWidth in $INVWIDTHPARAMS
      do
      dVal=${datasetDval##*,}
      dataset=${datasetDval%%,*}
      INDEX=$(($INDEX+1))
      MFILE=$NAME$INDEX.m
      SFILE=$NAME$INDEX.sh
      echo $MFILE
      echo $SFILE
      createMatlabFile "$MFILE" "'${dataset}'" "$dataNum" "$KERN" "$invWidth" "$dVal"
      createScriptFile "$SFILE" "$RUNNAME" "$MFILE"
    done
  done
done
echo Starting jobs ...
INDEX=$((0))
for dataNum in $DATANUM
  do
  for dataset in $DATASETS
    do
    for invWidth in $INVWIDTHPARAMS
      do
      INDEX=$(($INDEX+1))
      SFILE=$NAME$INDEX.sh
      qsub $SFILE
      sleep 0.1
    done
  done
done
echo Remember to clear the scripts later
