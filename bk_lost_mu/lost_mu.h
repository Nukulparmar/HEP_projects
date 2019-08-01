#ifndef lost_mu_H
#define lost_mu_H

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

class lost_mu : public NtupleVariables{

 public:
  lost_mu(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~lost_mu();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);
  bool     electron_match_photon(TLorentzVector);
  
  // Intialize histos here

  TH1D *ll_mu;
  TH1D *h_ht,*h_met,*h_lead_ph_pt,*h_njets,*h_el_size,*h_mu_size;
  TFile *oFile;
  
};
#endif

#ifdef lost_mu_cxx

void lost_mu::BookHistogram(const char *outFileName) {


  
  oFile = new TFile(outFileName, "recreate");
  // Define Histos here

  ll_mu = new TH1D("ll_mu","lost mu in b-jets and njets bins",6,1,7);
  h_ht = new TH1D("h_ht","HT after all the preselection",700,0,7000);
  h_met = new TH1D("h_met","MET after all the preselection",200,0,2000);
  h_lead_ph_pt = new TH1D("h_lead_ph_pt","Leading p_{T}^{#gamma} after all the preselection",200,0,2000);
  h_njets = new TH1D("h_njets","NJets after all the preselection",20,0,20);
  h_el_size = new TH1D("h_el_size","Number of Electron after all the preselection",5,0,5);
  h_mu_size = new TH1D("h_mu_size","Number of Muons after all the preselection",5,0,5);
  
}


lost_mu::lost_mu(const TString &inputFileList, const char *outFileName, const char* dataset) {
  
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

Bool_t lost_mu::FillChain(TChain *chain, const TString &inputFileList) {

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

Long64_t lost_mu::LoadTree(Long64_t entry) {
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

lost_mu::~lost_mu() { 
  
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  // From here the rest of the code is not executed and is giving the break segmentation error.
  oFile->Write();
  
  oFile->Close();

}

#endif

