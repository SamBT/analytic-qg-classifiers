#include <iostream>
#include <string.h>

#include "Pythia8/Pythia.h"
#include "EventHandler.h"

#include "boost/program_options.hpp"

using namespace Pythia8;
using std::atoi;
using std::atof;
using std::cout;
using std::endl;
namespace po = boost::program_options;

int main( int argc, char* argv[] ) {

  int nev;
  double CF, CA;
  int kernelOrder;
  string inFile;
  vector<string> configs;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("cf",po::value<double>(&CF)->default_value(1.33333),"quark casimir")
    ("ca",po::value<double>(&CA)->default_value(3.0),"gluon casimir")
    ("kernel,k",po::value<int>(&kernelOrder)->default_value(1),"kernel order")
    ("input,i",po::value<string>(&inFile),"input LHE file")
    ("config,c",po::value< vector<string> >(&configs),"input config files")
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  cout << "Called with nev = " << nev << ", CF = " << CF << ", CA = " << CA << ", kernel order = " << kernelOrder << endl;

  // initialize Pythia
  Pythia pythia;
  EventHandler Analysis(0.4,20);

  //Pythia settings not covered in DIRE config file
  //Use DIRE parton shower
  pythia.readString("PartonShowers:model = 3");
  // no substructure in e+e- beams
  pythia.readString("PDF:lepton = off");
  // set the quark casimir value
  pythia.settings.parm("DireColorQCD:CF", CF);
  // set the gluon casimir value
  pythia.settings.parm("DireColorQCD:CA", CA);
  //set the kernel order
  pythia.settings.mode("DireTimes:kernelOrder",kernelOrder);
  //turn off hadronization
  pythia.readString("HadronLevel:all = off");
  //Messing with the splitting functions
  pythia.readString("DireTimes:doGeneralizedKernel = on");

  //Read in LHE file
  string lhe = "Beams:LHEF = ";
  lhe += inFile;
  pythia.readString(lhe);

  for (int j = 0; j < configs.size(); j++) {
    pythia.readFile(configs[j]);
  }

  pythia.init();

  //Have to read config again to get pythia to "remember" some settings after init
  //e.g. init will set alphaSorder = 2 no matter what you tell it before init
  for (int j = 0; j < configs.size(); j++) {
    pythia.readFile(configs[j]);
  }

  pythia.init();

  pythia.settings.writeFile("run_settings.txt",true);

  cout << "---------------------------STARTING EVENT GEN--------------------------- \n \n \n" << endl;

  // Loop to generate events
  Analysis.Begin();
  for (int iev = 1; iev++) {
    // print a status update every 100 events
    if (iev % 100 == 0) cout << "Showered " << iev << " events." << endl;

    bool good = Analysis.AnalyzeEvent(iev,pythia);
    if (!good) break;
  }
  Analysis.End();

  return 0;
}
