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

mkdir ${pythdir}/${proc}/kernel${kern}
mkdir ${pythdir}/${proc}/kernel${kern}/${jobname}
rm -rf ${pythdir}/${proc}/kernel${kern}/${jobname}/output_${j}
mkdir ${pythdir}/${proc}/kernel${kern}/${jobname}/output_${j}
cd ${pythdir}/${proc}/kernel${kern}/${jobname}/output_${j}

module load root

echo "in dir $PWD"

cp /global/homes/s/sambt/Jets/optimal-classifiers/event-generation/src/GenEvts .

config_dir=/global/homes/s/sambt/Jets/optimal-classifiers/event-generation/DIRE-config-files
./GenEvts ${config_dir}/${proc}.cmnd ${nev} ${CF} ${CA} ${kern}

mv test.root ${pythdir}/${proc}/kernel${kern}/${jobname}/${filename}_${j}.root
rm -rf ${pythdir}/${proc}/kernel${kern}/${jobname}/output_${j}

cd ;
date
