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
  lost_el(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data",const char *year="year");
  ~lost_el();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *,const char *);
  void     BookHistogram(const char *);
  bool     electron_match_photon(TLorentzVector);
  void     fill_hist(TH1D *,TH1D *,int,double,int,double);
  void     fill_hist_2D(TH2D*,int,int,double,double);
  void     SR_hist(TH1D *,int,double,int,double,const char *,int);
  
  // Intialize histos here

  TH1D *total1,*total2,*fail_accept1,*fail_accept2,*fail_id1,*fail_id2,*fail_iso1,*fail_iso2,*one_lep_cr1,*one_lep_cr2,*fake_photon1,*fake_photon2;
  TH1D *h_st,*h_met,*h_lead_ph_pt[2],*h_lead_ph_eta[2],*h_njets,*h_el_size[2],*h_mu_size[2],*h_el_pt[2],*h_el_eta[2],*h_mu_pt[2],*h_mu_eta[2];
  TH1D *gen_el_pt,*gen_el_size[5],*gen_el_eta,*gen_mu_pt,*gen_mu_eta,*gen_mu_size[5],*gen_tau_size;
  TH1D *gen_ph_pt,*gen_ph_size[5],*gen_ph_eta;
  TH1D *nel[6],*nmu[6];
  TH1D *mindr_gen_rec_el,*mindr_reco_el_ph[2],*mindr_gen_el_reco_ph[2];
  TH2D *jet_ptbins,*e_ptbins,*mindr2D_el_jet,*mindr2D_genel_jet;
  TH1D *jetpt_inbins[9],*ept_inbins[9],*mindr1D_el_jet[9],*mindr1D_genel_jet[9];
  TH1D *lost_event1,*lost_event2;
  TH1D *h_mindr_ph_lep,*h_mindr_ph_qu,*h_mindr_lep_ph,*h_mindr_qu_ph,*h_mindr_goodph_lep;
  TH1D *h_madminphotonDR;
  TH1D *cr_pred1,*cr_pred2;
  TH1D *srbins_lostel,*srbins_cr;
  TH2D *METNJ_B0_E0,*METNJ_B1_E0,*METNJ_B0_E1,*METNJ_B1_E1,*TF_B0,*TF_B1;
  TFile *oFile;
  
};
#endif

#ifdef lost_el_cxx

void lost_el::BookHistogram(const char *outFileName) {

  
  oFile = new TFile(outFileName, "recreate");
  // Define Histos here
  Double_t metbins[]={0,100,150,5000};
  Double_t njetsbins_b0[]={0,2,3,4,6,7,20};
  Double_t njetsbins_b1[]={0,2,4,6,7,20};
  
  METNJ_B0_E0 = new TH2D("METNJ_B0_E0","SR for Bjets=0 in met and njet bins",3,metbins,6,njetsbins_b0); 
  METNJ_B0_E1 = new TH2D("METNJ_B0_E1","CR for Bjets=0 in met and njet bins",3,metbins,6,njetsbins_b0);
  METNJ_B1_E0 = new TH2D("METNJ_B1_E0","SR for Bjets>=1 in met and njet bins",3,metbins,5,njetsbins_b1);
  METNJ_B1_E1 = new TH2D("METNJ_B1_E1","CR for Bjets>=1 in met and njet bins",3,metbins,5,njetsbins_b1);
  TF_B0       = new TH2D("TF_B0","Transfer factors for BJets=0 in met and njet bins",3,metbins,5,njetsbins_b0);
  TF_B1       = new TH2D("TF_B1","Transfer factors for BJets>=1 in met and njet bins",3,metbins,5,njetsbins_b1);
  
  // srbins_cr = new TH1D("srbins_cr","1 e^- CR in SR_bins",25,1,26);
  // srbins_lostel = new TH1D("srbins_lostel","lost electron in SR_bins",25,1,26);
  
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
  fake_photon1 = new TH1D("fake_photon_1","fake photon in b-jets and njets bins",6,1,7);
  fake_photon2 = new TH1D("fake_photon_2","fake photon in all the bins",16,1,17);
  lost_event1 = new TH1D("lost_event_1","lost event in b-jets and njets bins",6,1,7);
  lost_event2 = new TH1D("lost_event_2","lost event in all the bins",16,1,17);
  cr_pred1 = new TH1D("cr_pred_1","1 lep cr in b-jets and njets bins",6,1,7);
  cr_pred2 = new TH1D("cr_pred_2","1 lep cr in all the bins",16,1,17);
  
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
  mindr_gen_rec_el = new TH1D("mindr_genel_reco_el","MinDr between gen and reco electron",500,0,5);
  mindr_gen_el_reco_ph[0] = new TH1D("mindr_gen_el_reco_ph_0","MinDr between gen e and reco #gamma",500,0,5);
  mindr_gen_el_reco_ph[1] = new TH1D("mindr_gen_el_reco_ph_1","MinDr between gen e and reco #gamma after the fake photon cut",500,0,5);
  mindr_reco_el_ph[0] = new TH1D("mindr_reco_el_ph_0","MinDr between reoc e and reco #gamma",500,0,5);
  mindr_reco_el_ph[1] = new TH1D("mindr_reco_el_ph_1","MinDr between reoc e and reco #gamma after the fake photon cut",500,0,5);
  jet_ptbins = new TH2D("jet_pt_bins","Jet Pts in the njet bins",9,1,10,200,0,600);
  e_ptbins = new TH2D("e_pt_bins","Electon Pts in the njet bins",9,1,10,200,0,500);
  mindr2D_el_jet = new TH2D("mindr2D_el_jet","mindr2D electrons and jets in njet bins",9,1,10,200,0,5);
  mindr2D_genel_jet = new TH2D("mindr2D_genel_jet","mindr2D gen electrons and jets in njet bins",9,1,10,200,0,5);
  char temp[200],temp2[200];
  for(int i=0;i<9;i++)
    { sprintf(temp,"jetpt_inbins_1D_%d",i+1);
      sprintf(temp2,"Jet Pt in the njet bin number %d",i+1);
      jetpt_inbins[i]= new TH1D(temp,temp2,200,0,400);
      sprintf(temp,"el_pt_inbins_1D_%d",i+1);
      sprintf(temp2,"Electrons Pt in the njet bin number %d",i+1);
      ept_inbins[i]= new TH1D(temp,temp2,200,0,400);
      sprintf(temp,"mindr1D_el_jet_%d",i+1);
      sprintf(temp2,"MinDr reco electron and jet in the njet bin number %d",i+1);
      mindr1D_el_jet[i]= new TH1D(temp,temp2,600,0,6);
      sprintf(temp,"mindr1D_genel_jet_%d",i+1);
      sprintf(temp2,"MinDr reco gen electron and jet in the njet bin number %d",i+1);
      mindr1D_genel_jet[i]= new TH1D(temp,temp2,600,0,6);
    }
  h_mindr_ph_lep = new TH1D("mindr_ph_lep","MinDr between the #gamma and lepton (e,#mu,#tau)",200,0,2);
  h_mindr_ph_qu  = new TH1D("mindr_ph_qu","MinDr between the #gamma and quarks (q,g)",200,0,2);
  h_mindr_lep_ph = new TH1D("mindr_lep_ph","MinDr between the lepton (e,#mu,#tau) and #gamma ",200,0,2);
  h_mindr_qu_ph  = new TH1D("mindr_qu_ph","MinDr between the quarks (q,g) and #gamma ",200,0,2);
  h_madminphotonDR = new TH1D("madphotonminDr","MinDr between the #gamma and quarks (q,g) using madMinPhotonDeltaR",200,0,2);
  h_mindr_goodph_lep = new TH1D("mindr_goodph_lep","MinDr between the good #gamma and lepton (e,#mu,#tau)",200,0,2);
}


lost_el::lost_el(const TString &inputFileList, const char *outFileName, const char* dataset, const char* year) {
  
  string nameData=dataset;//vvv
  string nameyear=year;
  TChain *tree = new TChain("PreSelection");
  //  tree = new TChain("TreeMaker2/PreSelection");//vvv
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << " and year " << year << std::endl;
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

