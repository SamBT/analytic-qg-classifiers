#!/bin/bash
name=$1
j=$2
kern=$3
cfg=$4


date
pythdir=/global/project/projectdirs/atlas/sambt/pythia83-output/optimal-classifiers

mkdir ${pythdir}/kernel${kern}
mkdir ${pythdir}/kernel${kern}/${name}
rm -rf ${pythdir}/kernel${kern}/${name}/output_${j}
mkdir ${pythdir}/kernel${kern}/${name}/output_${j}
cd ${pythdir}/kernel${kern}/${name}/output_${j}

module load root

echo "in dir $PWD"

cp /global/homes/s/sambt/Jets/optimal-classifiers/pythia83-evtgen/src/GenEvts .
source /global/homes/s/sambt/Jets/optimal-classifiers/pythia83-evtgen/src/setup.sh

./GenEvts $cfg

mv test.root ${pythdir}/${proc}/kernel${kern}/${name}/output_${j}.root
rm -rf ${pythdir}/${proc}/kernel${kern}/${name}/output_${j}

cd ;
date
