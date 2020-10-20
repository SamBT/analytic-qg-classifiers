#!/bin/bash
name=$1
j=$2
kern=$3
cfg=$4


date
pythdir=/global/project/projectdirs/atlas/sambt/Pythia_output/optimal-classifiers
mkdir ${pythdir}/${proc}

mkdir ${pythdir}/${proc}/kernel${kern}
mkdir ${pythdir}/${proc}/kernel${kern}/${name}
rm -rf ${pythdir}/${proc}/kernel${kern}/${name}/output_${j}
mkdir ${pythdir}/${proc}/kernel${kern}/${name}/output_${j}
cd ${pythdir}/${proc}/kernel${kern}/${name}/output_${j}

module load root

echo "in dir $PWD"

cp /global/homes/s/sambt/Jets/optimal-classifiers/event-generation/src/GenEvts .
source /global/homes/s/sambt/environment_scripts/setup_pythia_boost.sh

config_dir=/global/homes/s/sambt/Jets/optimal-classifiers/event-generation/DIRE-config-files
./GenEvts $cfg

mv test.root ${pythdir}/${proc}/kernel${kern}/${name}/output_${j}.root
rm -rf ${pythdir}/${proc}/kernel${kern}/${name}/output_${j}

cd ;
date
