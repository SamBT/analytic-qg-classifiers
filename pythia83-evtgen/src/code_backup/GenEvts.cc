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
  vector<string> configs;
  string outname;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("nev,n",po::value<int>(&nev)->default_value(1000),"number of events")
    ("cf",po::value<double>(&CF)->default_value(1.33333),"quark casimir")
    ("ca",po::value<double>(&CA)->default_value(3.0),"gluon casimir")
    ("kernel,k",po::value<int>(&kernelOrder)->default_value(1),"kernel order")
    ("config,c",po::value< vector<string> >(&configs),"config files")
    ("out,o",po::value<string>(&outname)->default_value("test.root"),"output root file name")
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

  //Read config file(s)
  for (int j = 0; j < configs.size(); j++) {
    cout << "reading config file : " << configs[j] << endl;
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
  Analysis.Begin(outname);
  for (int iev = 1; iev <= nev; iev++) {
    // print a status update every 100 events
    if (iev % 100 == 0) cout << "Generated " << iev << " events." << endl;

    Analysis.AnalyzeEvent(iev,pythia);
  }
  Analysis.End();

  return 0;
}
