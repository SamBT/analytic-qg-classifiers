//=============================================================================================
//
//   Pythia generation of e+e- > qqbar and e+e- > ggbar through an intermediate Higgs
//          based on a modified version of the main06.cc Pythia example script
//			and the dire05.cc Dire example script
//
//          Eric M. Metodiev, MIT, 2018 - 2019
//
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

  char config;
  int nev = 1000;
  double CF = 1.33333;
  double CA = 3.00000;
  if (argc < 2) {
    cout << "Too few arguments! Usage : ./events dire-config.cmnd nEvents CF CA" << endl;
    return 0;
  }
  else {
    //config = argv[1];
    CF = atof(argv[2]);
    CA = atof(argv[3]);
  }

  // initialize Pythia and Dire
  Pythia *pythia = new Pythia();
  Dire *dire = new Dire();
  EventHandler Analysis(0.4,20);

  //Pythia settings not covered in DIRE config file
  //pythia.readString("Main:numberOfEvents = "+argv[2])
  // no substructure in e+e- beams
  pythia->readString("PDF:lepton = off");
  // set the quark casimir value
  pythia->settings.parm("DireColorQCD:CF", CF);
  // set the gluon casimir value
  pythia->settings.parm("DireColorQCD:CA", CA);

  // initialize Dire
  dire->init(pythia, argv[1]);

  //pythia.init(); necessary?

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
