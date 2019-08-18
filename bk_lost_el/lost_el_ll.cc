#define lost_el_cxx
#include "lost_el.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>
#include <TMath.h>


using namespace std;
Int_t id=0;

int main(int argc, char* argv[])
{

  if (argc < 2) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  lost_el ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void lost_el::EventLoop(const char *data,const char *inputFileList)
{ 
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  TString s_data=data;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;

  double survived_events =0;

  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
    
      // ==============print number of events done == == == == == == == =
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
	cout << 10 * k << " %" <<endl;
      decade = k;
      // cout<<"j:"<<jentry<<" fcurrent:"<<fCurrent<<endl;
      // ===============read this entry == == == == == == == == == == == 
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
  
      Double_t wt = 35.9*1000*(Weight);

      vector<TLorentzVector> goodphotons;

      for(int i=0;i<Photons->size();i++)
	{ if(((*Photons)[i].Pt()>100)&&(TMath::Abs((*Photons)[i].Eta())<2.4)&&((*Photons_fullID)[i]==1)&&((*Photons_hasPixelSeed)[i]<0.000001))
	    { goodphotons.push_back((*Photons)[i]);
	      
	    }
	  
	}
      TLorentzVector goodphoton;
      sortTLorVec(&goodphotons);
      if(goodphotons.size()!=0)
	{goodphoton = goodphotons[0];}

      if(goodphoton.Pt()<100) continue; // remove events with 0 photons or pt<100
      if(MET<100) continue;
      
      
      vector<TLorentzVector> goodjets;
      double mindr = 100.0,st=0;
      int mindrindex = -1, gamma_matching_jet_index = -1;
      for(int i=0;i<Jets->size();i++)
	{ 
	  if((*Jets)[i].Pt()>30 && abs((*Jets)[i].Eta())<=2.5)
	    {
	      double dr = goodphoton.DeltaR((*Jets)[i]);
	      if(dr<mindr)
		{ mindr = dr;
		  mindrindex=i;
		}
	    }
	}
      for(int i =0;i<Jets->size();i++)
	{ if((*Jets)[i].Pt()>30 && abs((*Jets)[i].Eta())<=2.5)
	    {
	      if(!(mindr<0.3 && i==mindrindex))
		{ goodjets.push_back((*Jets)[i]);
		  st+=(*Jets)[i].Pt();
		}
	    }
	}
      if(mindr<0.3)
	{
	  gamma_matching_jet_index = mindrindex;
	  st+=goodphoton.Pt();
	}
      sortTLorVec(&goodjets);
      double dphi1=5.0,dphi2=5.0;
      if(goodjets.size()>0) dphi1 = abs(DeltaPhi(METPhi,goodjets[0].Phi()));
      if(goodjets.size()>1) dphi2 = abs(DeltaPhi(METPhi,goodjets[1].Phi()));
      if(gamma_matching_jet_index>=0 && ((*Jets)[gamma_matching_jet_index].Pt())/(goodphoton.Pt())<1.0) continue;
      if(gamma_matching_jet_index<0) continue;
      if(dphi1>0.3 && dphi2>0.3) continue;
      if((st<800 || goodphoton.Pt()<100) && (st<500 || goodphoton.Pt()<190)) continue;

    
      // Preselection cuts applied except the veto leptons and veto isotracks

      // Filling Some histos for keeping a check on the cuts

      h_lead_ph_pt[0]->Fill(goodphoton.Pt(),wt);
      h_lead_ph_eta[0]->Fill(goodphoton.Eta(),wt);
      h_met->Fill(MET,wt);
      h_st->Fill(st,wt);
      h_njets->Fill(goodjets.size(),wt);
      h_el_size[0]->Fill(Electrons->size(),wt);
      for(int i=0;i<Electrons->size();i++)
	{ h_el_pt[0]->Fill((*Electrons)[i].Pt(),wt);
	  h_el_eta[0]->Fill((*Electrons)[i].Eta(),wt);
	}
      
      h_mu_size[0]->Fill(Muons->size(),wt);
      for(int i=0;i<Muons->size();i++)
	{ h_mu_pt[0]->Fill((*Muons)[i].Pt(),wt);
	  h_mu_eta[0]->Fill((*Muons)[i].Eta(),wt);
	}
      
      


      
      // Getting the Gen Information

      // ----------------------  Gen level -------------------------
      vector<TLorentzVector> gen_el, gen_mu,gen_tau,gen_tau_el,gen_tau_mu;
     
      for(int i=0;i<GenParticles->size();i++)
	{ if((*GenParticles)[i].Pt()!=0)
	    {
	      if((abs((*GenParticles_PdgId)[i])==11) && (abs((*GenParticles_ParentId)[i])<=24) && ((*GenParticles_Status)[i]==1))
		{ 
		  gen_el.push_back((*GenParticles)[i]);
		  if(abs((*GenParticles_ParentId)[i])<=15)
		    {gen_tau_el.push_back((*GenParticles)[i]);}
		}
	      if((abs((*GenParticles_PdgId)[i])==13) && (abs((*GenParticles_ParentId)[i])<=24) && ((*GenParticles_Status)[i]==1))
		{ 
		  gen_mu.push_back((*GenParticles)[i]);
		  if(abs((*GenParticles_ParentId)[i])<=15)
		    {gen_tau_mu.push_back((*GenParticles)[i]);}
		}
	      if((abs((*GenParticles_PdgId)[i])==15) && (abs((*GenParticles_ParentId)[i])<=24) && ((*GenParticles_Status)[i]==1))
		{ 
		  gen_tau.push_back((*GenParticles)[i]);
		}
	      
	    }
	}
      
      if(gen_el.size()==0 && gen_mu.size()==0 && gen_tau.size()==0) continue; // rejecting hadronic decays

      h_lead_ph_pt[1]->Fill(goodphoton.Pt(),wt);
      h_lead_ph_eta[1]->Fill(goodphoton.Eta(),wt);
      
      h_el_size[1]->Fill(Electrons->size(),wt);
      for(int i=0;i<Electrons->size();i++)
	{ h_el_pt[1]->Fill((*Electrons)[i].Pt(),wt);
	  h_el_eta[1]->Fill((*Electrons)[i].Eta(),wt);
	}
      
      h_mu_size[1]->Fill(Muons->size(),wt);
      for(int i=0;i<Muons->size();i++)
	{ h_mu_pt[1]->Fill((*Muons)[i].Pt(),wt);
	  h_mu_eta[1]->Fill((*Muons)[i].Eta(),wt);
	}

      gen_el_size->Fill(gen_el.size(),wt);
      for(int i=0;i<gen_el.size();i++)
	{ gen_el_pt->Fill(gen_el[i].Pt(),wt);
	  gen_el_eta->Fill(gen_el[i].Eta(),wt);
	}
      
      gen_mu_size->Fill(gen_mu.size(),wt);
      for(int i=0;i<gen_mu.size();i++)
	{ gen_mu_pt->Fill(gen_mu[i].Pt(),wt);
	  gen_mu_eta->Fill(gen_mu[i].Eta(),wt);
	}
      
      
      
      fill_hist(*total1,*total2,BTags,MET,goodjets.size(),wt);
      
      
      
      
    }// loop over entries

  cout<<"Events survived = "<<survived_events<<endl;
}
	
  
bool lost_el::electron_match_photon(TLorentzVector photon)
{ for(int i=0;i<Electrons->size();i++)
    { if(photon.DeltaR((*Electrons)[i])<0.2)
	{ return true;
	}
    
      else
	{ return false;}

    }
}

void fill_hist(TH1D h1, TH1D h2,int btags,double met,int njets,double wt)
{
  h1.Fill("NJets_{=0}^{2-4}",0);
  h1.Fill("NJets_{#geq 1}^{2-4}",0);
  h1.Fill("NJets_{=0}^{5-6}",0);
  h1.Fill("NJets_{#geq 1}^{5-6}",0);
  h1.Fill("NJets_{=0}^{#geq 7}",0);
  h1.Fill("NJets_{#geq 1}^{#geq 7}",0);
  
  //int njets = goodjets.size();
  
  h2.Fill("NJets_{0}^{=2} & 100<MET<150",0);
  h2.Fill("NJets_{0}^{=2} & MET#geq 150",0);
  h2.Fill("NJets_{0}^{=3} & 100<MET<150",0);
  h2.Fill("NJets_{0}^{=3} & MET#geq 150",0);
  h2.Fill("NJets_{0}^{=4} & 100<MET<150",0);
  h2.Fill("NJets_{0}^{=4} & MET#geq 150",0);
  h2.Fill("NJets_{0}^{5-6} & 100<MET<150",0);
  h2.Fill("NJets_{0}^{5-6} & MET#geq 150",0);
  h2.Fill("NJets_{0}^{#geq7} & 100<MET<150",0);
  h2.Fill("NJets_{0}^{#geq7} & MET#geq 150",0);
  h2.Fill("NJets_{#geq 1}^{2-4} & 100<MET<150",0);
  h2.Fill("NJets_{#geq 1}^{2-4} & MET#geq 150",0);
  h2.Fill("NJets_{#geq 1}^{5-6} & 100<MET<150",0);
  h2.Fill("NJets_{#geq 1}^{5-6} & MET#geq 150",0);
  h2.Fill("NJets_{#geq 1}^{#geq 7} & 100<MET<150",0);
  h2.Fill("NJets_{#geq 1}^{#geq 7} & MET#geq 150",0);
  
  if(btags==0)
    { if(njets>=2 && njets<=4)
	{		  
	  h1.Fill("NJets_{=0}^{2-4}",wt);
	  if(njets==2)
	    {
	      if(met>=100 && met<150)
		{ h2.Fill("NJets_{0}^{=2} & 100<MET<150",wt);}
	      if(met>=150)
		{ h2.Fill("NJets_{0}^{=2} & MET#geq 150",wt);}
	    }
	  if(njets==3)
	    {
	      if(met>=100 && met<150)
		{ h2.Fill("NJets_{0}^{=3} & 100<MET<150",wt);}
	      if(met>=150)
		{ h2.Fill("NJets_{0}^{=3} & MET#geq 150",wt);}
	    }
	  if(njets==4)
	    {
	      if(met>=100 && met<150)
		{ h2.Fill("NJets_{0}^{=4} & 100<MET<150",wt);}
	      if(met>=150)
		{ h2.Fill("NJets_{0}^{=4} & MET#geq 150",wt);}
	    }
	}
      if(njets>=5 && njets<=6)
	{
	  h1.Fill("NJets_{=0}^{5-6}",wt);
	  if(met>=100 && met<150)
	    { h2.Fill("NJets_{0}^{5-6} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2.Fill("NJets_{0}^{5-6} & MET#geq 150",wt);}
	}
      if(njets>=7)
	{
	  h1.Fill("NJets_{=0}^{#geq 7}",wt);
	  if(met>=100 && met<150)
	    { h2.Fill("NJets_{0}^{#geq7} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2.Fill("NJets_{0}^{#geq7} & MET#geq 150",wt);}
	}

    }
  if(btags>=1)
    { if(njets>=2 && njets<=4)
	{
	  h1.Fill("NJets_{#geq 1}^{2-4}",wt);
	  if(met>=100 && met<150)
	    { h2.Fill("NJets_{#geq 1}^{2-4} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2.Fill("NJets_{#geq 1}^{2-4} & MET#geq 150",wt);}
	}
      if(njets>=5 && njets<=6)
	{
	  h1.Fill("NJets_{#geq 1}^{5-6}",wt);
	  if(met>=100 && met<150)
	    { h2.Fill("NJets_{#geq 1}^{5-6} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2.Fill("NJets_{#geq 1}^{5-6} & MET#geq 150",wt);}
	}
      if(njets>=7)
	{
	  h1.Fill("NJets_{#geq 1}^{#geq 7}",wt);
	  if(met>=100 && met<150)
	    { h2.Fill("NJets_{#geq 1}^{#geq 7} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2.Fill("NJets_{#geq 1}^{#geq 7} & MET#geq 150",wt);}
	}
	      
    }
}
