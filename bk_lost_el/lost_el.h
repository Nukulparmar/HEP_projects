#ifndef lost_el_H
#define lost_el_H

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

class lost_el : public NtupleVariables{

 public:
  lost_el(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~lost_el();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);
  bool     electron_match_photon(TLorentzVector);
  
  // Intialize histos here

  TH1D *ll,*ll_2;
  TH1D *h_ht,*h_met,*h_lead_ph_pt,*h_njets,*h_el_size,*h_mu_size,*h_dr_el_ph_1,*h_dr_el_ph_2,*h_dr_el_ph_3;
  TH1D *h_dr_ph_gen_el;
  TFile *oFile;
  
};
#endif

#ifdef lost_el_cxx

void lost_el::BookHistogram(const char *outFileName) {


  
  oFile = new TFile(outFileName, "recreate");
  // Define Histos here

  ll = new TH1D("ll_electron","lost electron in b-jets and njets bins",6,1,7);
  ll_2 = new TH1D("ll_electron_all","lost electron in all the bins",16,1,17);
  h_ht = new TH1D("h_ht","HT after all the preselection",700,0,7000);
  h_met = new TH1D("h_met","MET after all the preselection",200,0,2000);
  h_lead_ph_pt = new TH1D("h_lead_ph_pt","Leading p_{T}^{#gamma} after all the preselection",200,0,2000);
  h_njets = new TH1D("h_njets","NJets after all the preselection",20,0,20);
  h_el_size = new TH1D("h_el_size","Number of Electron after all the preselection",5,0,5);
  h_mu_size = new TH1D("h_mu_size","Number of Muons after all the preselection",5,0,5);
  h_dr_el_ph_1 = new TH1D("*h_dr_el_ph_1","dr between el and ph before any matching",100,0,5);
  h_dr_el_ph_2 = new TH1D("*h_dr_el_ph_2","dr between el and ph after reco el and reco ph matching",100,0,5);
  h_dr_el_ph_3 = new TH1D("*h_dr_el_ph_3","dr between el and ph after preselection cuts",100,0,5);
  h_dr_ph_gen_el = new TH1D("*h_dr_ph_gen_el","dr between gen el and ph after veto cuts",100,0,5);
  
}


lost_el::lost_el(const TString &inputFileList, const char *outFileName, const char* dataset) {
  
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

Bool_t lost_el::FillChain(TChain *chain, const TString &inputFileList) {

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

Long64_t lost_el::LoadTree(Long64_t entry) {
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

lost_el::~lost_el() { 
  
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  // From here the rest of the code is not executed and is giving the break segmentation error.
  oFile->Write();
  
  oFile->Close();

}

#endif

