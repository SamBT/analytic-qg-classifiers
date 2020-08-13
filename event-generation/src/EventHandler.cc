#include "EventHandler.h"

//Constructor
EventHandler::EventHandler(double radius = 0.4, double min_pt = 20) {
  outName = "test.root";
  pTmin = min_pt;
  jet_radius = radius;
  m_jet_def = new fastjet::JetDefinition(fastjet::antikt_algorithm, radius);
  m_jet_def_part = new fastjet::JetDefinition(fastjet::antikt_algorithm, radius);
}

//Destructor
EventHandler::~EventHandler() {
  delete m_jet_def;
  delete m_jet_def_part;
}

//Begin analysis, set up output file and TTree branches
void EventHandler::Begin() {
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

  //Loop over particles in event and store the final state particles and final (pre-hadronization) state partons
  for (int i = 0; i < pyth.event.size(); i++) {
    Pythia8::Particle part = pyth.event[i];
    int id = part.id();
    int idAbs = part.idAbs();

    //Skip non final-state particles and/or neutrinos
    if (!part.isFinal() && !part.isFinalPartonLevel()) continue;
    if (idAbs == 12 || idAbs == 14 || idAbs == 16) continue;

    fastjet::PseudoJet p(part.px(), part.py(), part.pz(), part.e());
    p.set_user_index(id);

    if (part.isFinal()) {
      particlesForJets.push_back(p);
    }
    else {
      partonsForJets.push_back(p);
    }
  } //End particle loop

  //Finding regular + parton jets
  fastjet::ClusterSequence cs(particlesForJets,*m_jet_def);
  fastjet::ClusterSequence cs_part(partonsForJets,*m_jet_def_part);
  vector<fastjet::PseudoJet> Jets = fastjet::sorted_by_pt(cs.inclusive_jets(pTmin));
  vector<fastjet::PseudoJet> partonJets = fastjet::sorted_by_pt(cs_part.inclusive_jets(pTmin));

  //Filling jet variables
  for (int i = 0; i < Jets.size(); i++) {
    if (nJetsFilled > max_njets) {cout << "Warning: Exceeded "<< max_njets << " jets!" << endl; continue;}
    njets++;
    jet_pt[nJetsFilled] = Jets[i].pt();
    jet_eta[nJetsFilled] = Jets[i].eta();
    jet_phi[nJetsFilled] = Jets[i].phi();
    jet_m[nJetsFilled] = Jets[i].m();
    jet_mult[nJetsFilled] = Jets[i].constituents().size();
    int type = JetType(Jets[i],partonsForJets,false);
    if (type <= 6) {
      nqjets++;
      is_qjet[nJetsFilled] = true;
      is_gjet[nJetsFilled] = false;
      qjet_pt[nQJetsFilled] = Jets[i].pt();
      qjet_eta[nQJetsFilled] = Jets[i].eta();
      qjet_phi[nQJetsFilled] = Jets[i].phi();
      qjet_m[nQJetsFilled] = Jets[i].m();
      qjet_mult[nQJetsFilled] = Jets[i].constituents().size();
      nQJetsFilled++;
    }
    else if (type == 21) {
      ngjets++;
      is_qjet[nJetsFilled] = false;
      is_gjet[nJetsFilled] = true;
      gjet_pt[nGJetsFilled] = Jets[i].pt();
      gjet_eta[nGJetsFilled] = Jets[i].eta();
      gjet_phi[nGJetsFilled] = Jets[i].phi();
      gjet_m[nGJetsFilled] = Jets[i].m();
      gjet_mult[nGJetsFilled] = Jets[i].constituents().size();
      nGJetsFilled++;
    }
    else {
      is_qjet[nJetsFilled] = false;
      is_gjet[nJetsFilled] = false;
    }
    nJetsFilled++;
  }

  //Filling parton-jet variables
  for (int i = 0; i < partonJets.size(); i++) {
    if (npJetsFilled > max_njets) {cout << "Warning: Exceeded "<< max_njets << " parton jets!" << endl; continue;}
    npjets++;
    pjet_pt[npJetsFilled] = partonJets[i].pt();
    pjet_eta[npJetsFilled] = partonJets[i].eta();
    pjet_phi[npJetsFilled] = partonJets[i].phi();
    pjet_m[npJetsFilled] = partonJets[i].m();
    pjet_mult[npJetsFilled] = partonJets[i].constituents().size();
    int type = JetType(partonJets[i],partonsForJets,true);
    if (type <= 6) {
      npqjets++;
      is_pqjet[npJetsFilled] = true;
      is_pgjet[npJetsFilled] = false;
      pqjet_pt[npQJetsFilled] = partonJets[i].pt();
      pqjet_eta[npQJetsFilled] = partonJets[i].eta();
      pqjet_phi[npQJetsFilled] = partonJets[i].phi();
      pqjet_m[npQJetsFilled] = partonJets[i].m();
      pqjet_mult[npQJetsFilled] = partonJets[i].constituents().size();
      npQJetsFilled++;
    }
    else if (type == 21) {
      npgjets++;
      is_pqjet[npJetsFilled] = false;
      is_pgjet[npJetsFilled] = true;
      pgjet_pt[npGJetsFilled] = partonJets[i].pt();
      pgjet_eta[npGJetsFilled] = partonJets[i].eta();
      pgjet_phi[npGJetsFilled] = partonJets[i].phi();
      pgjet_m[npGJetsFilled] = partonJets[i].m();
      pgjet_mult[npGJetsFilled] = partonJets[i].constituents().size();
      npGJetsFilled++;
    }
    else {
      is_pqjet[npJetsFilled] = false;
      is_pgjet[npJetsFilled] = false;
    }
    npJetsFilled++;
  }

  //Filling tree and finising up
  T->Fill();
  return;
}

int EventHandler::JetType(fastjet::PseudoJet jet, vector<fastjet::PseudoJet> partons, bool isPartonJet) {
  //Return |PDG ID| of highest-energy parton inside jet cone; determines if the jet is quark or gluon
  if (isPartonJet) {
    vector<fastjet::PseudoJet> ordered_partons = fastjet::sorted_by_E(jet.constituents());
    return abs(ordered_partons[0].user_index());
  }
  else {
    double maxE = 0;
    int id = -999;
    for (int i = 0; i < partons.size(); i++) {
      if (partons[i].delta_R(jet) < jet_radius && partons[i].e() > maxE) {
        id = abs(partons[i].user_index());
        maxE = partons[i].e();
      }
    }
    return id;
  }
}

void EventHandler::DeclareBranches() {
  T->Branch("EventNumber",&EventNumber,"EventNumber/I");
  T->Branch("pTmin",&pTmin,"pTmin/D");
  T->Branch("jetRadius",&jet_radius,"jet_radius/D");

  //Regular jets
  T->Branch("nJets",&njets,"njets/I");
  T->Branch("nQJets",&nqjets,"nqjets/I");
  T->Branch("nGJets",&ngjets,"ngjets/I");

  T->Branch("jet_pt",&jet_pt,"jet_pt[njets]/D");
  T->Branch("jet_eta",&jet_eta,"jet_eta[njets]/D");
  T->Branch("jet_phi",&jet_phi,"jet_phi[njets]/D");
  T->Branch("jet_m",&jet_m,"jet_m[njets]/D");
  T->Branch("jet_mult",&jet_mult,"jet_mult[njets]/I");
  T->Branch("is_qjet",&is_qjet,"is_qjet[njets]");
  T->Branch("is_gjet",&is_gjet,"is_gjet[njets]");

  T->Branch("qjet_pt",&qjet_pt,"qjet_pt[nqjets]/D");
  T->Branch("qjet_eta",&qjet_eta,"qjet_eta[nqjets]/D");
  T->Branch("qjet_phi",&qjet_phi,"qjet_phi[nqjets]/D");
  T->Branch("qjet_m",&qjet_m,"qjet_m[nqjets]/D");
  T->Branch("qjet_mult",&qjet_mult,"qjet_mult[nqjets]/I");

  T->Branch("gjet_pt",&gjet_pt,"gjet_pt[ngjets]/D");
  T->Branch("gjet_eta",&gjet_eta,"gjet_eta[ngjets]/D");
  T->Branch("gjet_phi",&gjet_phi,"gjet_phi[ngjets]/D");
  T->Branch("gjet_m",&gjet_m,"gjet_m[ngjets]/D");
  T->Branch("gjet_mult",&gjet_mult,"gjet_mult[ngjets]/I");

  //Parton jets
  T->Branch("npJets",&npjets,"npjets/I");
  T->Branch("npQJets",&npqjets,"npqjets/I");
  T->Branch("npGJets",&npgjets,"npgjets/I");

  T->Branch("pjet_pt",&pjet_pt,"pjet_pt[npjets]/D");
  T->Branch("pjet_eta",&pjet_eta,"pjet_eta[npjets]/D");
  T->Branch("pjet_phi",&pjet_phi,"pjet_phi[npjets]/D");
  T->Branch("pjet_m",&pjet_m,"pjet_m[npjets]/D");
  T->Branch("pjet_mult",&pjet_mult,"pjet_mult[npjets]/I");
  T->Branch("is_pqjet",&is_pqjet,"is_pqjet[npjets]");
  T->Branch("is_pgjet",&is_pgjet,"is_pgjet[npjets]");

  T->Branch("pqjet_pt",&pqjet_pt,"pqjet_pt[npqjets]/D");
  T->Branch("pqjet_eta",&pqjet_eta,"pqjet_eta[npqjets]/D");
  T->Branch("pqjet_phi",&pqjet_phi,"pqjet_phi[npqjets]/D");
  T->Branch("pqjet_m",&pqjet_m,"pqjet_m[npqjets]/D");
  T->Branch("pqjet_mult",&pqjet_mult,"pqjet_mult[npqjets]/I");

  T->Branch("pgjet_pt",&pgjet_pt,"pgjet_pt[npgjets]/D");
  T->Branch("pgjet_eta",&pgjet_eta,"pgjet_eta[npgjets]/D");
  T->Branch("pgjet_phi",&pgjet_phi,"pgjet_phi[npgjets]/D");
  T->Branch("pgjet_m",&pgjet_m,"pgjet_m[npgjets]/D");
  T->Branch("pgjet_mult",&pgjet_mult,"pgjet_mult[npgjets]/I");

  T->GetListOfBranches()->ls();

  return;
}

void EventHandler::ResetBranches() {
  EventNumber = -999;

  nJetsFilled = 0;
  nQJetsFilled = 0;
  nGJetsFilled = 0;
  njets = 0;
  nqjets = 0;
  ngjets = 0;

  npJetsFilled = 0;
  npQJetsFilled = 0;
  npGJetsFilled = 0;
  npjets = 0;
  npqjets = 0;
  npgjets = 0;

  for (int i = 0; i < max_njets; i++) {
    jet_pt[i] = -999;
    jet_eta[i] = -999;
    jet_phi[i] = -999;
    jet_m[i] = -999;
    jet_mult[i] = -999;
    is_qjet[i] = -999;
    is_gjet[i] = -999;

    qjet_pt[i] = -999;
    qjet_eta[i] = -999;
    qjet_phi[i] = -999;
    qjet_m[i] = -999;
    qjet_mult[i] = -999;

    gjet_pt[i] = -999;
    gjet_eta[i] = -999;
    gjet_phi[i] = -999;
    gjet_m[i] = -999;
    gjet_mult[i] = -999;

    pjet_pt[i] = -999;
    pjet_eta[i] = -999;
    pjet_phi[i] = -999;
    pjet_m[i] = -999;
    pjet_mult[i] = -999;
    is_pqjet[i] = -999;
    is_pgjet[i] = -999;

    pqjet_pt[i] = -999;
    pqjet_eta[i] = -999;
    pqjet_phi[i] = -999;
    pqjet_m[i] = -999;

    pgjet_pt[i] = -999;
    pgjet_eta[i] = -999;
    pgjet_phi[i] = -999;
    pgjet_m[i] = -999;
  }

  return;
}
