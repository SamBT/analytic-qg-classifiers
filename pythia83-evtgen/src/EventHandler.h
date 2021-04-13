#ifndef EventHandler_H
#define EventHandler_H

#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <set>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "Pythia8/Pythia.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

using namespace std;

class EventHandler {
  private:
    //Tree & ROOT output objects
    TTree *T;
    TFile *F;
    string outName;

    //Event variables
    int EventNumber;

    //Jet variables
    static const int max_njets = 20;
    double pTmin = 20;
    double jet_radius;

    int nJetsFilled;
    int nQJetsFilled;
    int nGJetsFilled;

    int njets;
    double jet_pt[max_njets];
    double jet_eta[max_njets];
    double jet_phi[max_njets];
    double jet_m[max_njets];
    int jet_mult[max_njets];

    double lead_constit_pt[50];
    double lead_constit_eta[50];
    double lead_constit_phi[50];
    double lead_constit_e[50];
    int lead_constit_id[50];

    //Fastjet objects
    fastjet::JetDefinition *m_jet_def; //Regular anti-kT jets

  public:
    EventHandler (double radius, double min_pt);
    ~EventHandler ();
    void Begin(string out);
    void AnalyzeEvent(int iEvt, Pythia8::Pythia& pyth);
    void End();
    void DeclareBranches();
    void ResetBranches();
    void SetOutName(string fname){
      outName = fname;
    }
};

#endif
