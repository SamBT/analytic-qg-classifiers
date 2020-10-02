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
    bool is_qjet[max_njets];
    bool is_gjet[max_njets];

    //Quark jet variables
    int nqjets;
    double qjet_pt[max_njets];
    double qjet_eta[max_njets];
    double qjet_phi[max_njets];
    double qjet_m[max_njets];
    int qjet_mult[max_njets];

    //Gluon jet variables
    int ngjets;
    double gjet_pt[max_njets];
    double gjet_eta[max_njets];
    double gjet_phi[max_njets];
    double gjet_m[max_njets];
    int gjet_mult[max_njets];

    //Parton Jet Variables
    int npJetsFilled;
    int npQJetsFilled;
    int npGJetsFilled;

    int npjets;
    double pjet_pt[max_njets];
    double pjet_eta[max_njets];
    double pjet_phi[max_njets];
    double pjet_m[max_njets];
    int pjet_mult[max_njets];
    double plead_constit_pt[100];
    double plead_constit_eta[100];
    double plead_constit_phi[100];
    double plead_constit_e[100];
    int plead_constit_id[100];
    bool is_pqjet[max_njets];
    bool is_pgjet[max_njets];

    //Quark parton-jet variables
    int npqjets;
    double pqjet_pt[max_njets];
    double pqjet_eta[max_njets];
    double pqjet_phi[max_njets];
    double pqjet_m[max_njets];
    int pqjet_mult[max_njets];
    double pqlead_constit_pt[100];
    double pqlead_constit_eta[100];
    double pqlead_constit_phi[100];
    double pqlead_constit_e[100];

    //Gluon parton-jet variables
    int npgjets;
    double pgjet_pt[max_njets];
    double pgjet_eta[max_njets];
    double pgjet_phi[max_njets];
    double pgjet_m[max_njets];
    int pgjet_mult[max_njets];
    double pglead_constit_pt[100];
    double pglead_constit_eta[100];
    double pglead_constit_phi[100];
    double pglead_constit_e[100];

    //Fastjet objects
    fastjet::JetDefinition *m_jet_def; //Regular anti-kT jets
    fastjet::JetDefinition *m_jet_def_part; //Parton-jets

  public:
    EventHandler (double radius, double min_pt);
    ~EventHandler ();
    void Begin();
    void AnalyzeEvent(int iEvt, Pythia8::Pythia& pyth);
    int JetType(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> partons);
    void End();
    void DeclareBranches();
    void ResetBranches();
    void SetOutName(string fname){
      outName = fname;
    }
};

#endif
