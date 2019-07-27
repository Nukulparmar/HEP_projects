#ifndef AnalyzeLightBSM_H
#define AnalyzeLightBSM_H

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

class AnalyzeLightBSM : public NtupleVariables{

 public:
  AnalyzeLightBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~AnalyzeLightBSM();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);

  // Intialize histos here

  TH1D *h_gen_photon_pt,*h_gen_lead_photon_pt,*h_lead_photon_pt,*h_photon_pt,*h_photon_eta,*h_photon_phi,*h_HT,*h_Njets,*h_MET,*h_Nbjets,*h_MHT,*h_isoeltracks,*h_isomutracks,*h_isopitracks,*h_dr_photon_jets,*h_dr_el_fake_photon,*h_el_pt,*h_mu_pt,*h_eff,*h_dr_jet_fake_el,*h_dr_jet_fake_ph;

  TH2D *genel_vs_recphoton,*eff_vs_jetpt,*el_ph_fake_corr,*genjet_vs_recel,*jet_el_fake_corr,*genjet_vs_recph,*jet_ph_fake_corr;
  
  
  TFile *oFile;
  
};
#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName) {


  
  oFile = new TFile(outFileName, "recreate");
  // Define Histos here

  h_gen_photon_pt = new TH1D("h_gen_photon_pt","gen #gamma pt",200,0,2000);
  h_gen_lead_photon_pt = new TH1D("h_gen_lead_photon_pt","gen leading #gamma pt",200,0,1400);
  h_lead_photon_pt = new TH1D("h_lead_photon_pt","leading #photon pt",200,0,1600);
  h_photon_pt = new TH1D("h_photon_pt","#photon pt",200,0,1600);
  h_photon_eta = new TH1D("h_photon_eta","#photon eta",100,-3,3);
  h_photon_phi = new TH1D("h_photon_phi","#photon phi",200,-4,4);
  h_HT = new TH1D("h_HT","HT without removing the overlap with #gamma",100,0,5000);
  h_Njets = new TH1D("h_Njets","NJets without removing the overlap with #gamma",20,0,20); 
  h_MET = new TH1D("h_MET","MET",100,0,1000);
  h_Nbjets = new TH1D("h_Nbjets","Number of jets with Btag",20,0,20);
  h_MHT = new TH1D("h_MHT","MHT",100,0,800);
  h_isoeltracks = new TH1D("iso_el_tracks","Number of isolated electron tracks",5,0,5);
  h_isomutracks = new TH1D("iso_mu_tracks","Number of isolated muon tracks",5,0,5);
  h_isopitracks = new TH1D("iso_pi_tracks","Number of isolated pion tracks",5,0,5);
  h_dr_photon_jets = new TH1D("dr_photon_jets","DeltaR(photon,jets)",100,0,10);
  h_dr_el_fake_photon = new TH1D("dr_el_fake_photon","DeltaR(genElectron,recoPhoton) for Electron faking photon",100,0,6);
  genel_vs_recphoton = new TH2D("genel_vs_recphoton","gen Electron faking as reco photons",100,0,2000,100,0,2000);
  h_el_pt = new TH1D("h_el_pt","electron pt",200,0,2000);
  h_mu_pt = new TH1D("h_mu_pt","muon pt",200,0,2000);
  h_eff = new TH1D("eff","Mismeasurement eff (GenJet-Jet)/GenJet",100,0,1);
  eff_vs_jetpt = new TH2D("eff_vs_jetpt","Mismeasurement eff vs jet pt",200,0,800,100,0,1);
  el_ph_fake_corr = new TH2D("el_ph_fake_corr","Electron Faking Photon correlation",100,0,800,40,0,4);
  h_dr_jet_fake_el = new TH1D("dr_jet_fake_el","DeltaR(genJet,electron) for Jets faking electrons",100,0,6);
  genjet_vs_recel = new TH2D("genjet_vs_recel","genjet faking as rec e pt correlation",100,0,500,100,0,500);
  jet_el_fake_corr =new TH2D("jet_el_fake_corr","Jet faking Electron correlation",100,0,600,40,0,2);
  h_dr_jet_fake_ph = new TH1D("dr_jet_fake_ph","DeltaR(genJet,photon) for Jets faking electrons",100,0,6);
  genjet_vs_recph = new TH2D("genjet_vs_recph","genjet faking as rec #gamma pt correlation",100,0,600,100,0,500);
  jet_ph_fake_corr =new TH2D("jet_ph_fake_corr","Jet faking Photon correlation",100,0,600,40,0,2);
}


AnalyzeLightBSM::AnalyzeLightBSM(const TString &inputFileList, const char *outFileName, const char* dataset) {
  
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

Bool_t AnalyzeLightBSM::FillChain(TChain *chain, const TString &inputFileList) {

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

Long64_t AnalyzeLightBSM::LoadTree(Long64_t entry) {
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

AnalyzeLightBSM::~AnalyzeLightBSM() { 
  
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  // From here the rest of the code is not executed and is giving the break segmentation error.
  oFile->Write();
  
  oFile->Close();

}

#endif

