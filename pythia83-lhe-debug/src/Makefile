CXX=g++
BOOST=/global/project/projectdirs/atlas/sambt/MC_Generators/boost-install
BOOSTINC=/global/project/projectdirs/atlas/sambt/MC_Generators/boost-install/include
BOOSTLIB=/global/project/projectdirs/atlas/sambt/MC_Generators/boost-install/lib
DIRE_CONFIG=/global/project/projectdirs/atlas/sambt/MC_Generators/DIRE-2.004/bin/dire-config

GenEvts : GenEvts.cc EventHandler.cc
	$(CXX) -o $@ $^ \
	-I/global/project/projectdirs/atlas/sambt/MC_Generators/pythia8303//include -L/global/project/projectdirs/atlas/sambt/MC_Generators/pythia8303//lib -Wl,-rpath /global/project/projectdirs/atlas/sambt/MC_Generators/pythia8303//lib -lpythia8 -I/global/project/projectdirs/atlas/sambt/MC_Generators/fastjet-install//include -L/global/project/projectdirs/atlas/sambt/MC_Generators/fastjet-install//lib -Wl,-rpath /global/project/projectdirs/atlas/sambt/MC_Generators/fastjet-install//lib -lfastjet -DFASTJET3 -I/global/project/projectdirs/atlas/sambt/MC_Generators/hepmc2.06.11-install//include -L/global/project/projectdirs/atlas/sambt/MC_Generators/hepmc2.06.11-install//lib -Wl,-rpath /global/project/projectdirs/atlas/sambt/MC_Generators/hepmc2.06.11-install//lib -lHepMC -DHEPMC2 -I/global/project/projectdirs/atlas/sambt/MC_Generators/LHAPDF-6.3.0-install//include -L/global/project/projectdirs/atlas/sambt/MC_Generators/LHAPDF-6.3.0-install//lib -Wl,-rpath /global/project/projectdirs/atlas/sambt/MC_Generators/LHAPDF-6.3.0-install//lib -lLHAPDF -DLHAPDF6 -ldl -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -std=c++11 \
	`root-config --cflags --glibs` \
	-I$(BOOSTINC) \
	-L$(BOOSTLIB) -lboost_program_options
