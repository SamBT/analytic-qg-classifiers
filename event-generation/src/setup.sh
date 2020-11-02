setup_PYTHIA() {
    export PYTHIA8LOCATION=/global/project/projectdirs/atlas/sambt/MC_Generators/pythia8244/
    export PYTHIA8DATA=${PYTHIA8LOCATION}/share/Pythia8/xmldoc
    export LD_LIBRARY_PATH=${PYTHIA8LOCATION}/lib/:$LD_LIBRARY_PATH
}

setup_fastjet() {
    export FASTJETLOCATION=/global/project/projectdirs/atlas/sambt/MC_Generators/fastjet-install
    export LD_LIBRARY_PATH=/global/project/projectdirs/atlas/sambt/MC_Generators/fastjet-install/lib:$LD_LIBRARY_PATH
}

setup_boost() {
    export BOOSTINCDIR=/global/project/projectdirs/atlas/sambt/MC_Generators/boost-install/include
    export BOOSTLIBLOCATION=/global/project/projectdirs/atlas/sambt/MC_Generators/boost-install/lib
    export LD_LIBRARY_PATH=${BOOSTLIBLOCATION}:$LD_LIBRARY_PATH
}

setup_PYTHIA
setup_fastjet
setup_boost
