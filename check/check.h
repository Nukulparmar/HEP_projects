#ifndef check_H
#define check_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2.h"
#include <TProfile.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TGraph.h"

class check : public NtupleVariables{

 public:
  check(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~check();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);
  bool     electron_match_photon(TLorentzVector);
  double   mindr(TLorentzVector,vector<TLorentzVector>);
  int      boolvalue(bool);
  //  void     fill_histo(TH1D,TH1D,vector<TLorentzVector>,double);
  
  // Intialize histos here

  TH1D *num_events,*btags;
  TH1D *el_mediumId,*el_miniIso,*el_passIso,*el_size,*el_isotrack,*el_nelectrons,*el_tightId;
  TH1D *mu_mediumId,*mu_miniIso,*mu_passIso,*mu_size,*mu_isotrack,*mu_nmuons,*mu_tightId;
  TH1D *ph_elfakes,*ph_pixelseed,*ph_isEB,*ph_passelveto,*ph_pfgammaiso,*ph_fullId;
  TH1D *pion_isotrack;
  
  TH1D *lead_ph_pt,*lead_ph_eta,*all_ph_pt,*all_ph_eta,*all_el_pt,*all_el_eta,*all_mu_pt,*all_mu_eta,*gen_el_pt,*gen_el_eta,*gen_mu_pt,*gen_mu_eta,*gen_tau_pt,*gen_tau_eta,*genel_size1,*genel_size2,*genmu_size1,*genmu_size2,*gentau_size1,*gentau_size2,*gentop_size,*gentop_pt,*gentop_eta;

  TH1D *genel_fromtau,*genmu_fromtau,*gentau_had;
  TH1D *mindr_recoel_recoph,*mindr_recoel_genel,*mindr_recoph_genel,*mindr_recoph_genph;
  TH2D *genel_recoel_eta,*genel_recoph_eta,*genph_recoph_eta,*recoel_recoph_eta;
  
  
  TFile *oFile;
  
};
#endif

#ifdef check_cxx

void check::BookHistogram(const char *outFileName) {


  
  oFile = new TFile(outFileName, "recreate");
  // Define Histos here

  num_events = new TH1D("numevents","Number of events",200,0,170000);
  btags = new TH1D("BTags","BTags",6,0,6);
  lead_ph_pt = new TH1D("lead_ph_pt","lead #gamma PT",200,0,2000);
  lead_ph_eta = new TH1D("lead_ph_eta","lead #gamma Eta",100,0,5);
  all_ph_pt = new TH1D("all_ph_pt","All reco #gamma PT",200,0,2000);
  all_ph_eta = new TH1D("all_ph_eta","All reco #gamma Eta",100,0,5);
  ph_elfakes = new TH1D("ph_electronFakes","Photon -> Electron Fakes",2,0,2);
  ph_pixelseed = new TH1D("ph_pixelseed","Photon hass pixel seed",2,0,2);
  ph_isEB = new TH1D("ph_isEB","Photon is in EB",2,0,2);
  ph_passelveto = new TH1D("ph_passelectronveto","Photon passing electron veto",2,0,2);
  ph_pfgammaiso= new TH1D("ph_pfGammaIso","Photons Pf Gamma Iso",200,0,10);
  ph_fullId = new TH1D("ph_fullId","Photon full Id",2,0,2);
  

  all_el_pt = new TH1D("all_el_pt","All reco e^{-} PT",200,0,1000);
  all_el_eta = new TH1D("all_el_eta","All reco e^{-} Eta",100,0,5);
  el_size = new TH1D("el_size","Size of reco e{-}",6,0,6);
  el_mediumId = new TH1D("Electrons_mediumID","Electrons Medium ID",3,0,3);
  el_miniIso = new TH1D("Electrons_miniIso","Electrons MiniIso",200,0,20);
  el_passIso = new TH1D("Electrons_passIso","Electrons Pass Iso",2,0,2);
  el_tightId = new TH1D("Electrons_tightId","Electrons tight Id",2,0,2);
  el_isotrack = new TH1D("IsoElectronTracks","Isolated Electron tracks",5,0,5);
  el_nelectrons = new TH1D("NElectrons","Electrons after Isolation cut (Check this)",5,0,5);
  
  
  all_mu_pt = new TH1D("all_mu_pt","All reco #mu PT",200,0,1000);
  all_mu_eta = new TH1D("all_mu_eta","All reco #mu Eta",100,0,5);
  mu_size = new TH1D("Muons_size","Size of reco #mu",6,0,6);
  mu_mediumId = new TH1D("Muons_mediumID","Muons Medium ID",3,0,3);
  mu_miniIso = new TH1D("Muons_miniIso","Muons MiniIso",200,0,20);
  mu_passIso = new TH1D("Muons_passIso","Muons Pass Iso",2,0,2);
  mu_tightId = new TH1D("Muons_tightId","Muons tight Id",2,0,2);
  mu_isotrack = new TH1D("IsoMuonsTracks","Isolated Muons tracks",5,0,5);
  mu_nmuons = new TH1D("NMuons","Muons after Isolation cut (Check this)",5,0,5);
  
  
  /* all_tau_pt = new TH1D("all_tau_pt","All reco #tau PT",200,0,1000); */
  /* all_tau_eta = new TH1D("all_tau_eta","All reco #tau Eta",100,0,5); */


  pion_isotrack= new TH1D("IsoPionTracks","Isolated Pion Tracks",5,0,5);
  
  gen_el_pt = new TH1D("gen_el_pt","All gen e^{-} PT",200,0,1000);
  gen_el_eta = new TH1D("gen_el_eta","All gen e^{-} Eta",100,0,5);
  genel_size1 = new TH1D("genel_size1","Gen electron size",6,0,6);
  genel_size2 = new TH1D("genel_size2","Gen electron size from GenElectrons branch",6,0,6);
  genel_fromtau = new TH1D("genel_fromtau","gen electron form #tau",2,0,2);
  
  gen_mu_pt = new TH1D("gen_mu_pt","All gen #mu PT",200,0,1000);
  gen_mu_eta = new TH1D("gen_mu_eta","All gen #mu Eta",100,0,5);
  genmu_size1 = new TH1D("genmu_size1","Gen #mu size",6,0,6);
  genmu_size2 = new TH1D("genmu_size2","Gen #mu size from GenMuons branch",6,0,6);
  genmu_fromtau = new TH1D("genmu_fromtau","gen #mu form #tau",2,0,2);
  
  gen_tau_pt = new TH1D("gen_tau_pt","All gen #tau PT",200,0,1000);
  gen_tau_eta = new TH1D("gen_tau_eta","All gen #tau Eta",100,0,5);
  gentau_size1 = new TH1D("gentau_size1","Gen #tau size",6,0,6);
  gentau_size2 = new TH1D("gentau_size2","Gen #tau size from GenTaus branch",6,0,6);
  gentau_had = new TH1D("gentau_had","Gen #tau hadronic",2,0,2);

  gentop_pt = new TH1D("gentop_PT","Gen Top Pt",200,0,5000);
  gentop_eta = new TH1D("gentop_eta","Gen Top eta",100,0,5);
  gentop_size = new TH1D("gentop_size","Gen Top size",3,0,3);
  mindr_recoel_recoph = new TH1D("mindr_reco_el_reco_ph","MinDR between Reco e and Reco #gamam",200,0,6);
  mindr_recoel_genel = new TH1D("mindr_recoel_genel","MinDR between Reco e and Gen e",200,0,6);
  mindr_recoph_genel = new TH1D("mindr_recoph_genel","MinDR between Reco #gamma and Gen e",200,0,6);
  mindr_recoph_genph = new TH1D("mindr_recoph_genph","MinDR between Reco #gamma and Gen #gamma",200,0,6);

  genel_recoel_eta = new TH2D("genel_recoel_eta","Gen e vs reco e Eta",200,0,5,200,0,5);
  genel_recoph_eta = new TH2D("genel_recoph_eta","Gen e vs reco #gamma Eta",200,0,5,200,0,5);
  genph_recoph_eta = new TH2D("genph_recoph_eta","Gen #gamma vs reco #gamma Eta",200,0,5,200,0,5);
  recoel_recoph_eta = new TH2D("recoph_recoph_eta","Gen e vs reco #gamma Eta",200,0,5,200,0,5);
}


check::check(const TString &inputFileList, const char *outFileName, const char* dataset) {
  
  string nameData=dataset;//vvv
  TChain *tree = new TChain("PreSelection");
  //  tree = new TChain("TreeMaker2/PreSelection");//vvv
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }
  
  /* if(nameData!="signalH") nameData="BG"; */
  /* if(nameData=="signalH") nameData="signal"; */
  //cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  NtupleVariables::Init(tree,nameData);
  
  BookHistogram(outFileName);
  
}

Bool_t check::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    
    return kFALSE;
  }
  
  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
}

Long64_t check::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    
    //Notify();
  }
  return centry;
}

check::~check() { 
  
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  // From here the rest of the code is not executed and is giving the break segmentation error.
  oFile->Write();
  
  oFile->Close();

}

#endif

