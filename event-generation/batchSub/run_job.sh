#!/bin/bash
jobname=$1
filename=$2
j=$3
proc=$4
nev=$5
CF=$6
CA=$7
kern=$8


date
pythdir=/global/project/projectdirs/atlas/sambt/Pythia_output/optimal-classifiers
mkdir ${pythdir}/${proc}

rm -rf ${pythdir}/${proc}/${jobname}_${j}
mkdir ${pythdir}/${proc}/${jobname}_${j}
cd ${pythdir}/${proc}/${jobname}_${j}

module load root

cp /global/homes/s/sambt/Jets/optimal-classifiers/event-generation/src/GenEvts .

config_dir=/global/homes/s/sambt/Jets/optimal-classifiers/event-generation/DIRE-config-files
./GenEvts ${config_dir}/${proc}.cmnd ${nev} ${CF} ${CA} ${kern}

mv test.root ${filename}_${j}.root

cd ;
date
