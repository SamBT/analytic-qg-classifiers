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

using namespace Pythia8;

int main( int argc, char* argv[] ) {

  // initialize Pythia and Dire
  Pythia pythia;
  Dire dire;
  dire.initSettings(pythia);

  // parse the command line inputs or use the defaults
  int nev = 10000;
  int process = 0;
  int ecm = 200;
  double CF = 1.33333;
  double CA = 3.00000;
  for (int i = 0; i < argc; i++){
    if (strcmp(argv[i],"-nev") == 0) nev = atoi(argv[i+1]);
    if (strcmp(argv[i],"-proc")== 0) process = atoi(argv[i+1]);
    if (strcmp(argv[i],"-ecm") == 0) ecm = atoi(argv[i+1]);
    if (strcmp(argv[i],"-CF")  == 0) CF = atof(argv[i+1]);
    if (strcmp(argv[i],"-CA")  == 0) CA = atof(argv[i+1]);
  }

  // output filename
  char filename[100];
  strcpy(filename, "../events/");
  
  // turn off quark mass
  pythia.readString("1:m0=000");
  
  switch (process) {
    // e+ e- > q qbar through an intermediate Higgs boson
    case 0:
      pythia.readString("HiggsSM:ffbar2H = on");
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfAny = 1");
      strcat(filename, "qq_events");
      break;

    // e+ e- > g g through an intermediate Higgs boson
    case 1:
      pythia.readString("HiggsSM:ffbar2H = on");
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfAny = 21");
      strcat(filename, "gg_events");
      break;

    default:
      cerr << "Invalid process number!";
      return 1;
  }

  // set the quark casimir value
  pythia.settings.parm("DireColorQCD:CF", CF);
  strcat(filename, "_CF");
  std::stringstream ss1;
  ss1 << std::fixed << std::setprecision(0) << 100.0*CF;
  strcat(filename, ss1.str().c_str());
  
  // set the gluon casimir value
  pythia.settings.parm("DireColorQCD:CA", CA);
  strcat(filename, "_CA");
  std::stringstream ss2;
  ss2 << std::fixed << std::setprecision(0) << 100.0*CA;
  strcat(filename, ss2.str().c_str());

  // initialize Dire
  dire.init(pythia, argv[1]);

  // set the output file
  strcat(filename, "_ko0.txt");
  ofstream outstream(filename);
  std::cout << filename << endl;
  
  // no substructure in e+e- beams
  pythia.readString("PDF:lepton = off");

  // e+e- collisions at a specified center of mass
  pythia.readString("Beams:idA =  11");
  pythia.readString("Beams:idB = -11");
  pythia.settings.parm("Beams:eCM", ecm);
  pythia.init();
  
  // Loop to generate events
  for (int iEvent = 0; iEvent < nev;) {
    if (!pythia.next()) continue;

    // iterate over the particles in the event
    for (int i = 0; i < pythia.event.size(); i++) {
      
      const Pythia8::Particle & particle = pythia.event[i];
      const int idabs = particle.idAbs();
      
      // final state only, no neutrinos
      if (!particle.isFinal() || idabs == 12 || idabs == 14 || idabs == 16) continue;
      
      // write out the particle four-momentum
      outstream << particle.e() <<" "<< particle.px() <<" "<< particle.py() <<" "<< particle.pz() << " " << particle.m() <<" "<< particle.id();
      outstream << endl;
    }  
    outstream << endl;
    iEvent += 1;
    
    // print a status update every 100 jets
    if (iEvent % 100 == 0) 
      cout << "Generated " << iEvent << " events." << endl;
  }
  return 0;
}
