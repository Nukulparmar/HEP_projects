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
  void     fill_hist(TH1D *,TH1D *,int,double,int,double);

  // Intialize histos here

  TH1D *total1,*total2,*fail_accept1,*fail_accept2,*fail_id1,*fail_id2,*fail_iso1,*fail_iso2,*one_lep_cr1,*one_lep_cr2,*fake_photon1,*fake_photon2,*hadtau1,*hadtau2;
  TH1D *h_st,*h_met,*h_lead_ph_pt[2],*h_lead_ph_eta[2],*h_njets,*h_el_size[2],*h_mu_size[2],*h_el_pt[2],*h_el_eta[2],*h_mu_pt[2],*h_mu_eta[2];
  TH1D *gen_el_pt,*gen_el_size[5],*gen_el_eta,*gen_mu_pt,*gen_mu_eta,*gen_mu_size[5],*gen_tau_size;
  TH1D *gen_ph_pt,*gen_ph_size[5],*gen_ph_eta;
  TH1D *nel[6],*nmu[6];
  TH1D *mindr_gen_rec_mu,*mindr_reco_el_ph[2],*mindr_gen_el_reco_ph[2];
 
  TFile *oFile;
  
  
};
#endif

#ifdef lost_mu_cxx

void lost_mu::BookHistogram(const char *outFileName) {


  
  oFile = new TFile(outFileName, "recreate");
  // Define Histos here

 
  total1 = new TH1D("total1","lost electron in b-jets and njets bins",6,1,7);
  total2 = new TH1D("total2","lost electron in all the bins",16,1,17);
  fail_accept1 = new TH1D("fail_accept_1","Fail Acceptance in b-jets and njets bins",6,1,7);
  fail_accept2 = new TH1D("fail_accept_2","Fail Acceptance in all the bins",16,1,17);
  fail_id1 = new TH1D("fail_id_1","Fail Id in b-jets and njets bins",6,1,7);
  fail_id2 = new TH1D("fail_id_2","Fail Id in all the bins",16,1,17);
  fail_iso1 = new TH1D("fail_iso_1","Fail Iso in b-jets and njets bins",6,1,7);
  fail_iso2 = new TH1D("fail_iso_2","Fail Iso in all the bins",16,1,17);
  one_lep_cr1 = new TH1D("one_lep_cr_1","1 lep cr in b-jets and njets bins",6,1,7);
  one_lep_cr2 = new TH1D("one_lep_cr_2","1 lep cr in all the bins",16,1,17);
  fake_photon1 = new TH1D("fake_photon_1","1 lep cr in b-jets and njets bins",6,1,7);
  fake_photon2 = new TH1D("fake_photon_2","1 lep cr in all the bins",16,1,17);
  hadtau1 = new TH1D("hadtau1","Hadronic #tau in b-jets and njets bins",6,1,7);
  hadtau2 = new TH1D("hadtau2","Hadronic  #tau in all the bins",16,1,17);
  
  
  h_st = new TH1D("h_ht","HT after all the preselection",700,0,7000);
  h_met = new TH1D("h_met","MET after all the preselection",200,0,2000);
  h_lead_ph_pt[0] = new TH1D("h_lead_ph_pt0","Leading p_{T}^{#gamma} after all the preselection",200,0,2000);
  h_lead_ph_eta[0] = new TH1D("h_lead_ph_eta0","Leading #gamma #eta after all the preselection",200,-6,6);
  h_lead_ph_pt[1] = new TH1D("h_lead_ph_p1t","Leading p_{T}^{#gamma} after veto hadronic",200,0,2000);
  h_lead_ph_eta[1] = new TH1D("h_lead_ph_eta1","Leading #gamma #eta after veto hadronic",200,-6,6);
  h_njets = new TH1D("h_njets","NJets after all the preselection",20,0,20);
  
  h_el_size[0] = new TH1D("h_el_size0","Number of Electron after all the preselection",5,0,5);
  h_el_pt[0] = new TH1D("h_el_pt0","Reco e pt after preselection",200,0,1600);
  h_el_eta[0] = new TH1D("h_el_eta0","Reco e #eta after preselection",200,-6,6);
  h_el_size[1] = new TH1D("h_el_size1","Number of Electron after veto hadronic",5,0,5);
  h_el_pt[1] = new TH1D("h_el_pt1","Reco e pt after veto hadronic",200,0,1600);
  h_el_eta[1] = new TH1D("h_el_eta1","Reco e #eta after veto hadronic",200,-6,6);

  nel[0] = new TH1D("nel_0","NElectron after preselection",5,0,5);
  nel[1] = new TH1D("nel_1","NElectron after veto gen hadronic",5,0,5);
  nel[2] = new TH1D("nel_2","NElectron after fail acceptance",5,0,5);
  nel[3] = new TH1D("nel_3","NElectron after fail Id",5,0,5);
  nel[4] = new TH1D("nel_4","NElectron after fail Iso",5,0,5);
  nel[5] = new TH1D("nel_5","NElectron after all",5,0,5);

  nmu[0] = new TH1D("nmu_0","Nmuon after preselection",5,0,5);
  nmu[1] = new TH1D("nmu_1","Nmuon after veto gen hadronic",5,0,5);
  nmu[2] = new TH1D("nmu_2","Nmuon after fail acceptance",5,0,5);
  nmu[3] = new TH1D("nmu_3","Nmuon after fail Id",5,0,5);
  nmu[4] = new TH1D("nmu_4","Nmuon after fail Iso",5,0,5);
  nmu[5] = new TH1D("nmu_5","Nmuon after all",5,0,5);
  
  h_mu_size[0] = new TH1D("h_mu_size0","Number of Muons after all the preselection",5,0,5);
  h_mu_pt[0] = new TH1D("h_mu_pt0","Reco #mu pt after preselection",200,0,1600);
  h_mu_eta[0] = new TH1D("h_mu_eta0","Reco #mu #eta after preselection",200,-6,6);
  h_mu_size[1] = new TH1D("h_mu_size1","Number of Muons after veto hadronic",5,0,5);
  h_mu_pt[1] = new TH1D("h_mu_pt1","Reco #mu pt after veto hadronic",200,0,1600);
  h_mu_eta[1] = new TH1D("h_mu_eta1","Reco #mu #eta after veto hadronic",200,-6,6);

  gen_el_pt = new TH1D("gen_el_pt","Gen el pt after preselection and veto hadronic",200,0,1600);
  gen_el_eta = new TH1D("gen_el_eta","Gen el #eta after preselection and veto hadronic",200,-6,6);
  gen_el_size[0] = new TH1D("gen_el_size_0","Gen el size after preslection and veto hadronic",5,0,5);
  gen_el_size[1] = new TH1D("gen_el_size_1","Gen el size after fail acceptance",5,0,5);
  gen_el_size[2] = new TH1D("gen_el_size_2","Gen el size after fail Id",5,0,5);
  gen_el_size[3] = new TH1D("gen_el_size_3","Gen el size after fail Iso",5,0,5);
  gen_el_size[4] = new TH1D("gen_el_size_4","Gen el size after all",5,0,5);

  
  gen_mu_pt = new TH1D("gen_mu_pt","Gen mu pt after preselection and veto hadronic",200,0,1600);
  gen_mu_eta = new TH1D("gen_mu_eta","Gen mu #eta after preselection and veto hadronic",200,-6,6);
  gen_mu_size[0] = new TH1D("gen_mu_size_0","Gen mu size after preslection and veto hadronic",5,0,5);
  gen_mu_size[1] = new TH1D("gen_mu_size_1","Gen mu size after fail acceptance",5,0,5);
  gen_mu_size[2] = new TH1D("gen_mu_size_2","Gen mu size after fail Id",5,0,5);
  gen_mu_size[3] = new TH1D("gen_mu_size_3","Gen mu size after fail Iso",5,0,5);
  gen_mu_size[4] = new TH1D("gen_mu_size_4","Gen mu size after all",5,0,5);

  gen_ph_pt = new TH1D("gen_ph_pt","Gen ph pt after preselection and veto hadronic",200,0,1600);
  gen_ph_eta = new TH1D("gen_ph_eta","Gen ph #eta after preselection and veto hadronic",200,-6,6);
  gen_ph_size[0] = new TH1D("gen_ph_size_0","Gen ph size after preslection and veto hadronic",5,0,5);
  gen_ph_size[1] = new TH1D("gen_ph_size_1","Gen ph size after fail acceptance",5,0,5);
  gen_ph_size[2] = new TH1D("gen_ph_size_2","Gen ph size after fail Id",5,0,5);
  gen_ph_size[3] = new TH1D("gen_ph_size_3","Gen ph size after fail Iso",5,0,5);
  gen_ph_size[4] = new TH1D("gen_ph_size_4","Gen ph size after all",5,0,5);
  
  gen_tau_size = new TH1D("gen_tau_size","Gen tau size after preslection and veto hadronic",5,0,5);
  mindr_gen_rec_mu = new TH1D("mindr_genel_reco_el","MinDr between gen and reco electron",500,0,5);
  mindr_gen_el_reco_ph[0] = new TH1D("mindr_gen_el_reco_ph_0","MinDr between gen e and reco #gamma",500,0,5);
  mindr_gen_el_reco_ph[1] = new TH1D("mindr_gen_el_reco_ph_1","MinDr between gen e and reco #gamma after the fake photon cut",500,0,5);
  mindr_reco_el_ph[0] = new TH1D("mindr_reco_el_ph_0","MinDr between reoc e and reco #gamma",500,0,5);
  mindr_reco_el_ph[1] = new TH1D("mindr_reco_el_ph_1","MinDr between reoc e and reco #gamma after the fake photon cut",500,0,5);

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

