#!/bin/bash
#SBATCH -A m3246
#SBATCH -q shared
#SBATCH -t 00:30:00
#SBATCH -C haswell
#SBATCH -c 4

name=$1
j=$2
kern=$3
cfg=$4

export pythdir=/global/project/projectdirs/atlas/sambt/pythia83-output/optimal-classifiers

mkdir -p ${pythdir}/kernel${kern}
mkdir -p ${pythdir}/kernel${kern}/${name}
rm -rf ${pythdir}/kernel${kern}/${name}/output_${j}
mkdir -p ${pythdir}/kernel${kern}/${name}/output_${j}

module load root

cp /global/homes/s/sambt/Jets/optimal-classifiers/pythia83-evtgen/src/GenEvts ${pythdir}/kernel${kern}/${name}/output_${j}
source /global/homes/s/sambt/Jets/optimal-classifiers/pythia83-evtgen/src/setup.sh

tstart=`date +%s`

cd ${pythdir}/kernel${kern}/${name}/output_${j}

srun ./GenEvts $cfg

mv test.root ${pythdir}/${proc}/kernel${kern}/${name}/output_${j}.root
cd ${pythdir}/${proc}/kernel${kern}/${name}
rm -rf ${pythdir}/${proc}/kernel${kern}/${name}/output_${j}

tstop=`date +%s`
tot=`expr $((tstop - tstart)) / 60`
secs=`expr $((tstop - tstart)) % 60`
echo "Execution took $tot minutes $secs seconds"
