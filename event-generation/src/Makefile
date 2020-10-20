CXX=g++
BOOST=/global/project/projectdirs/atlas/sambt/MC_Generators/boost-install
BOOSTINC=/global/project/projectdirs/atlas/sambt/MC_Generators/boost-install/include
BOOSTLIB=/global/project/projectdirs/atlas/sambt/MC_Generators/boost-install/lib
DIRE_CONFIG=/global/project/projectdirs/atlas/sambt/MC_Generators/DIRE-2.004/bin/dire-config

GenEvts : GenEvts.cc EventHandler.cc
	$(CXX) -o $@ $^ \
	`$(DIRE_CONFIG) --all-libs` \
	`root-config --cflags --glibs` \
	-I$(BOOSTINC) \
	-L$(BOOSTLIB) -lboost_program_options
