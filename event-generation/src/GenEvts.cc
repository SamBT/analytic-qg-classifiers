//=============================================================================================
//
//   Pythia generation of e+e- > qqbar and e+e- > ggbar through an intermediate Higgs
//          based on a modified version of the main06.cc Pythia example script
//			and the dire05.cc Dire example script
//    g++ events.cc -o events `/Users/metodiev/bin/dire-config --all-libs`
//    g++ events.cc -o events `/Users/metodiev/pythia8235/examples/direforpythia-66862c4022ad6958140db46680af723723acfe26/trunk/bin/dire-config --all-libs`
//    ./events dire.cmnd -nev 1000 -proc 1 -ecm 200 -CF 1.333333 -CA 3.000000
//
//=============================================================================================

#include <iostream>
#include <string.h>

#include "Pythia8/Pythia.h"
#include "Dire/Dire.h"
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
  vector<string> configs;
  
  /*if (argc < 5) {
    cout << "Too few arguments! Usage : ./events dire-config.cmnd nEvents CF CA kernelOrder" << endl;
    return 0;
  }
  else {
    //config = argv[1];
    nev = atoi(argv[2]);
    CF = atof(argv[3]);
    CA = atof(argv[4]);
    kernelOrder = atoi(argv[5]);
    }*/

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("nev,n",po::value<int>(&nev)->default_value(1000),"number of events")
    ("cf",po::value<double>(&CF)->default_value(1.33333),"quark casimir")
    ("ca",po::value<double>(&CA)->default_value(3.0),"gluon casimir")
    ("kernel,k",po::value<int>(&kernelOrder)->default_value(1),"kernel order")
    ("config,c",po::value< vector<string> >(&configs),"config files")
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  cout << "Called with nev = " << nev << ", CF = " << CF << ", CA = " << CA << ", kernel order = " << kernelOrder << endl;

  // initialize Pythia and Dire
  Pythia pythia;
  Dire dire;
  EventHandler Analysis(0.4,20);
  dire.initSettings(pythia);

  //Pythia settings not covered in DIRE config file
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

  //Read config file(s)
  for (int j = 0; j < configs.size(); j++) {
    cout << "reading config file : " << configs[j] << endl;
    pythia.readFile(configs[j]);
  }

  dire.init(pythia);
  pythia.init();

  cout << "---------------------------STARTING EVENT GEN--------------------------- \n \n \n" << endl;

  // Loop to generate events
  Analysis.Begin();
  for (int iev = 1; iev <= nev; iev++) {
    // print a status update every 100 events
    if (iev % 100 == 0) cout << "Generated " << iev << " events." << endl;

    Analysis.AnalyzeEvent(iev,pythia);
  }
  Analysis.End();

  return 0;
}
