#include "EventHandler.h"

//Constructor
EventHandler::EventHandler(double radius = 0.4, double min_pt = 20){
  pTmin = min_pt;
  jet_radius = radius;
  m_jet_def = new fastjet::JetDefinition(fastjet::antikt_algorithm, radius);
}

//Destructor
EventHandler::~EventHandler() {
  delete m_jet_def;
}

//Begin analysis, set up output file and TTree branches
void EventHandler::Begin(string out) {
  outName = out;
  F = new TFile(outName.c_str(),"RECREATE");
  T = new TTree("EventTree","Generated Events");

  DeclareBranches();
  ResetBranches();

  return;
}

//Write output and end analysis
void EventHandler::End() {
  T->Write();
  F->Close();
  return;
}

void EventHandler::AnalyzeEvent(int iEvt, Pythia8::Pythia& pyth) {
  if (!pyth.next()) return;

  ResetBranches();

  EventNumber = iEvt;
  vector<fastjet::PseudoJet> particlesForJets; //particle list for final-state jets
  vector<fastjet::PseudoJet> partonsForJets; //parton list for "truth"/parton jets

  //Loop over particles in event and store the final state partons
  for (int i = 0; i < pyth.event.size(); i++) {
    Pythia8::Particle part = pyth.event[i];
    int id = part.id();
    int idAbs = part.idAbs();
    
    //Skip non final-state particles and/or neutrinos
    if (!part.isFinal() || idAbs == 12 || idAbs == 14 || idAbs == 16) continue;

    fastjet::PseudoJet p(part.px(), part.py(), part.pz(), part.e());
    p.set_user_index(id);

    //Keep final-state partons for jet clustering
    particlesForJets.push_back(p);
  } //End particle loop
  EventPartonMult = particlesForJets.size();

  //Finding regular + parton jets
  fastjet::ClusterSequence cs(particlesForJets,*m_jet_def);
  vector<fastjet::PseudoJet> Jets = fastjet::sorted_by_pt(cs.inclusive_jets(pTmin));

  //Filling jet variables
  for (int i = 0; i < Jets.size(); i++) {
    if (nJetsFilled > max_njets) {cout << "Warning: Exceeded "<< max_njets << " jets!" << endl; continue;}
    njets++;
    jet_pt[nJetsFilled] = Jets[i].pt();
    jet_eta[nJetsFilled] = Jets[i].eta();
    jet_phi[nJetsFilled] = Jets[i].phi();
    jet_m[nJetsFilled] = Jets[i].m();
    jet_mult[nJetsFilled] = Jets[i].constituents().size();
    if (i == 0) {
      for (int k = 0; k < Jets[i].constituents().size(); k++) {
        lead_constit_pt[k] = Jets[i].constituents()[k].pt();
        lead_constit_eta[k] = Jets[i].constituents()[k].eta();
        lead_constit_phi[k] = Jets[i].constituents()[k].phi();
        lead_constit_e[k] = Jets[i].constituents()[k].e();
        lead_constit_id[k] = Jets[i].constituents()[k].user_index();
      }
    }
    nJetsFilled++;
  }

  //Filling tree and finising up
  T->Fill();
  return;
}

void EventHandler::DeclareBranches() {
  T->Branch("EventNumber",&EventNumber,"EventNumber/I");
  T->Branch("pTmin",&pTmin,"pTmin/D");
  T->Branch("jetRadius",&jet_radius,"jet_radius/D");
  T->Branch("eventPartonMult",&EventPartonMult,"EventPartonMult/I");

  //Regular jets
  T->Branch("nJets",&njets,"njets/I");

  T->Branch("jet_pt",&jet_pt,"jet_pt[njets]/D");
  T->Branch("jet_eta",&jet_eta,"jet_eta[njets]/D");
  T->Branch("jet_phi",&jet_phi,"jet_phi[njets]/D");
  T->Branch("jet_m",&jet_m,"jet_m[njets]/D");
  T->Branch("jet_mult",&jet_mult,"jet_mult[njets]/I");
  
  T->Branch("lead_constit_pt",&lead_constit_pt,"plead_constit_pt[200]/D");
  T->Branch("lead_constit_eta",&lead_constit_eta,"plead_constit_eta[200]/D");
  T->Branch("lead_constit_phi",&lead_constit_phi,"plead_constit_phi[200]/D");
  T->Branch("lead_constit_e",&lead_constit_e,"plead_constit_e[200]/D");
  T->Branch("lead_constit_id",&lead_constit_id,"plead_constit_id[200]/I");

  T->GetListOfBranches()->ls();

  return;
}

void EventHandler::ResetBranches() {
  EventNumber = -999;

  nJetsFilled = 0;
  njets = 0;
  EventPartonMult = 0;

  for (int i = 0; i < max_njets; i++) {
    jet_pt[i] = -999;
    jet_eta[i] = -999;
    jet_phi[i] = -999;
    jet_m[i] = -999;
    jet_mult[i] = -999;
  }

  for (int i = 0; i < 200; i++) {
    lead_constit_pt[i] = 0;
    lead_constit_eta[i] = 0;
    lead_constit_phi[i] = 0;
    lead_constit_e[i] = 0;
    lead_constit_id[i] = 0;
  }

  return;
}
