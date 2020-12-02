#include <iostream>
#include <string.h>

#include "Pythia8/Pythia.h"

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

  // initialize Pythia
  Pythia pythia;
  //Pythia settings not covered in DIRE config file
  //Use DIRE parton shower
  pythia.readString("PartonLevel:all = off");
  pythia.readString("HadronLevel:all = off");
  // no substructure in e+e- beams
  pythia.readString("PDF:lepton = off");
  // set the quark casimir value
  pythia.settings.parm("DireColorQCD:CF", CF);
  // set the gluon casimir value
  pythia.settings.parm("DireColorQCD:CA", CA);
  //set the kernel order
  pythia.settings.mode("DireTimes:kernelOrder",kernelOrder);
  //turn off hadronization
  //Messing with the splitting functions
  pythia.readString("DireTimes:doGeneralizedKernel = on");

  //Create LHE writer object
  LHAupFromPYTHIA8 myLHA(&pythia.process, &pythia.info);
  myLHA.openLHEF("hardProcess.lhe");

  //Read config file(s)
  for (int j = 0; j < configs.size(); j++) {
    cout << "reading config file : " << configs[j] << endl;
    pythia.readFile(configs[j]);
  }

  pythia.init();
  // Store initialization info in the LHAup object.
  myLHA.setInit();
  // Write out this initialization info on the file.
  myLHA.initLHEF();
  //Have to read config again to get pythia to "remember" some settings after init
  //e.g. init will set alphaSorder = 2 no matter what you tell it before init

  for (int j = 0; j < configs.size(); j++) {
    pythia.readFile(configs[j]);
  }

  // Loop to generate events
  for (int iev = 1; iev <= nev; iev++) {
    // print a status update every 100 events
    if (iev % 100 == 0) cout << "Generated " << iev << " events." << endl;
    // Generate an event.
    pythia.next();
    // Store event info in the LHAup object.
    myLHA.setEvent();
    // Write out this event info on the file.
    // With optional argument (verbose =) false the file is smaller.
    myLHA.eventLHEF();
  }
  // Update the cross section info based on Monte Carlo integration during run.
  myLHA.updateSigma();
  // Write endtag. Overwrite initialization info with new cross sections.
  myLHA.closeLHEF(true);

  return 0;
}
