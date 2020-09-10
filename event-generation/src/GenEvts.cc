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

using namespace Pythia8;
using std::atoi;
using std::atof;
using std::cout;
using std::endl;

int main( int argc, char* argv[] ) {

  int nev = 1000;
  double CF = 1.33333;
  double CA = 3.00000;
  int kernelOrder = 3;
  if (argc < 5) {
    cout << "Too few arguments! Usage : ./events dire-config.cmnd nEvents CF CA kernelOrder" << endl;
    return 0;
  }
  else {
    //config = argv[1];
    nev = atoi(argv[2]);
    CF = atof(argv[3]);
    CA = atof(argv[4]);
    kernelOrder = atoi(argv[5]);
  }

  // initialize Pythia and Dire
  Pythia pythia;
  Dire dire;
  EventHandler Analysis(0.4,20);
  //dire.init(pythia, argv[1]);
  dire.initSettings(pythia);

  //Pythia settings not covered in DIRE config file
  // no substructure in e+e- beams
  pythia.readString("PDF:lepton = off");
  pythia.readString("HadronLevel:all = off"); // don't need hadronization
  // set the quark casimir value
  pythia.settings.parm("DireColorQCD:CF", CF);
  // set the gluon casimir value
  pythia.settings.parm("DireColorQCD:CA", CA);
  //set the kernel order
  pythia.settings.mode("DireTimes:kernelOrder",kernelOrder);

  dire.init(pythia, argv[1]);
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
