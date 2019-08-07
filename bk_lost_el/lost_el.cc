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

  int survived_events =0;

  
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
      for(int i=0;i<Electrons->size();i++)
	{	  h_dr_el_ph_1->Fill(goodphoton.DeltaR((*Electrons)[i]),wt);         }
	    
	    
      if(electron_match_photon(goodphoton)) continue; // don't want electron matched to good photon
      for(int i=0;i<Electrons->size();i++)
	{	  h_dr_el_ph_2->Fill(goodphoton.DeltaR((*Electrons)[i]),wt);         }

      
      if(Muons->size()>0) continue; // veto muons
      if(Electrons->size()>1) continue; // 1 or no(lost) reco electrons
      if(isoMuonTracks!=0 || isoPionTracks!=0) continue;

      double mt_ele=0,mt_pho=0,mt_elepho=0;
      if(Electrons->size()==1){
	mt_ele=sqrt(2*(*Electrons)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*Electrons)[0].Phi()))));
	if(mt_ele>100) continue;
	if( ((*Electrons)[0].Pt() < 10) || abs((*Electrons)[0].Eta()) > 2.5 ) continue;
      }
      mt_pho=sqrt(2*goodphoton.Pt()*MET*(1-cos(DeltaPhi(METPhi,goodphoton.Phi()))));
      if(Electrons->size()==1) mt_elepho=sqrt(2*(goodphoton+(*Electrons)[0]).Pt()*MET*(1-cos(DeltaPhi(METPhi,(goodphoton+(*Electrons)[0]).Phi()))));   // MHT cut defined in the paper for isolated tracks


      vector<TLorentzVector> goodjets;
      double mindr = 100.0,ht=0;
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
		  ht+=(*Jets)[i].Pt();
		}
	    }
	}
      if(mindr<0.3)
	{
	  gamma_matching_jet_index = mindrindex;
	  ht+=goodphoton.Pt();
	}

      sortTLorVec(&goodjets);

      double dphi1=5.0,dphi2=5.0;
      if(goodjets.size()>0) dphi1 = abs(DeltaPhi(METPhi,goodjets[0].Phi()));
      if(goodjets.size()>1) dphi2 = abs(DeltaPhi(METPhi,goodjets[1].Phi()));

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

      if(Electrons->size()==0)
	{
	  if(isoMuonTracks!=0 || isoElectronTracks!=0 || isoPionTracks!=0) continue;
	  if(gen_el.size()==0) continue;
	  // if(gen_mu.size()!=0) continue;
	}
      if(gen_mu.size()!=0) continue; // rejecting all events with gen muons
      sortTLorVec(&gen_el);

      // checking if the photon is real or fake

      bool realphoton = true;
      int match_el=0,match_p=0;
      double mindr_ph_genobj = 100,match_el_pt=0.0,match_ph_pt=0.0;
      for(int i=0;i<GenParticles->size();i++)
	{ if((*GenParticles)[i].Pt()!=0)
	    { double dr1=goodphoton.DeltaR((*GenParticles)[i]);
	      h_dr_ph_gen_el->Fill(dr1,wt);
	      if(dr1<0.2 && abs((*GenParticles_PdgId)[i])==11 && abs((*GenParticles_ParentId)[i])<=24)
		{ match_el=1;
		  match_el_pt=(*GenParticles)[i].Pt();
		}
	      if(mindr_ph_genobj > dr1) { mindr_ph_genobj = dr1;}
	    }
	}

      for(int i=0;i<GenParticles->size();i++)
	{ if((*GenParticles)[i].Pt()!=0)
	    { double dr1=goodphoton.DeltaR((*GenParticles)[i]);
	      if(dr1<0.2 && abs((*GenParticles_PdgId)[i])==22)
		{ match_p=1;
		  match_ph_pt=(*GenParticles)[i].Pt();
		}
	    }
	}

      if(Electrons->size()==0)
	{
	  if(match_el==0) realphoton = true;
	  else if(match_el==1 && match_p==0) realphoton=false;
	  else if(match_el==1 && match_p==1)
	    { if(abs(goodphoton.Pt() - match_ph_pt) > abs(goodphoton.Pt()-match_el_pt))
		{ realphoton=false;}
	    }
	  else realphoton=true;
	}
      if(!realphoton) continue;

      // electron matching with jets and lost in a jet ( no need to check this as reco electron collection is defined as isolated electrons

      double mindr_el_jet = 100.0;
      int index_jet_match_el = -1;

      if(Electrons->size()==1)
      	{ for(int i=0;i<Jets->size();i++)
      	    {
      	      if((*Jets)[i].Pt()>30 && abs((*Jets)[i].Eta())<=2.5)
      		{
      		  double dr = (*Electrons)[0].DeltaR((*Jets)[i]);
      		  if(dr<mindr_el_jet)
      		    { mindr_el_jet = dr;
      		      index_jet_match_el = i;}
      		}
      	    }

      	   if(index_jet_match_el>=0 && ((*Jets)[index_jet_match_el].Pt())/((*Electrons)[0].Pt())<1.0) continue;
      	   if(index_jet_match_el<0) continue;
      	}
      if(gamma_matching_jet_index>=0 && ((*Jets)[gamma_matching_jet_index].Pt())/(goodphoton.Pt())<1.0) continue;
      if(gamma_matching_jet_index<0) continue;


      if(MET>100 && goodjets.size()>=2 && (dphi1>0.3 && dphi2 >0.3) && ht>100 && goodphoton.Pt()>100)
	{ survived_events+=1;
	  h_ht->Fill(ht,wt);
	  h_met->Fill(MET,wt);
	  h_lead_ph_pt->Fill(goodphoton.Pt(),wt);
	  h_njets->Fill(goodjets.size(),wt);
	  h_el_size->Fill(Electrons->size(),wt);
	  h_mu_size->Fill(Muons->size(),wt);
	  for(int i=0;i<Electrons->size();i++)
	    {	  h_dr_el_ph_3->Fill(goodphoton.DeltaR((*Electrons)[i]),wt);         }


	  ll->Fill("NJets_{=0}^{2-4}",0);
	  ll->Fill("NJets_{#geq 1}^{2-4}",0);
	  ll->Fill("NJets_{=0}^{5-6}",0);
	  ll->Fill("NJets_{#geq 1}^{5-6}",0);
	  ll->Fill("NJets_{=0}^{#geq 7}",0);
	  ll->Fill("NJets_{#geq 1}^{#geq 7}",0);

	  int njets = goodjets.size();

	  ll_2->Fill("NJets_{0}^{=2} & 100<MET<150",0);
	  ll_2->Fill("NJets_{0}^{=2} & MET#geq 150",0);
	  ll_2->Fill("NJets_{0}^{=3} & 100<MET<150",0);
	  ll_2->Fill("NJets_{0}^{=3} & MET#geq 150",0);
	  ll_2->Fill("NJets_{0}^{=4} & 100<MET<150",0);
	  ll_2->Fill("NJets_{0}^{=4} & MET#geq 150",0);
	  ll_2->Fill("NJets_{0}^{5-6} & 100<MET<150",0);
	  ll_2->Fill("NJets_{0}^{5-6} & MET#geq 150",0);
	  ll_2->Fill("NJets_{0}^{#geq7} & 100<MET<150",0);
	  ll_2->Fill("NJets_{0}^{#geq7} & MET#geq 150",0);
	  ll_2->Fill("NJets_{#geq 1}^{2-4} & 100<MET<150",0);
	  ll_2->Fill("NJets_{#geq 1}^{2-4} & MET#geq 150",0);
	  ll_2->Fill("NJets_{#geq 1}^{5-6} & 100<MET<150",0);
	  ll_2->Fill("NJets_{#geq 1}^{5-6} & MET#geq 150",0);
	  ll_2->Fill("NJets_{#geq 1}^{#geq 7} & 100<MET<150",0);
	  ll_2->Fill("NJets_{#geq 1}^{#geq 7} & MET#geq 150",0);

	  if(BTags==0)
	    { if(goodjets.size()>=2 && goodjets.size()<=4)
		{		  
		  ll->Fill("NJets_{=0}^{2-4}",1);
		  if(njets==2)
		    {
		      if(MET>=100 && MET<150)
			{ ll_2->Fill("NJets_{0}^{=2} & 100<MET<150",1);}
		      if(MET>=150)
			{ ll_2->Fill("NJets_{0}^{=2} & MET#geq 150",1);}
		    }
		  if(njets==3)
		    {
		      if(MET>=100 && MET<150)
			{ ll_2->Fill("NJets_{0}^{=3} & 100<MET<150",1);}
		      if(MET>=150)
			{ ll_2->Fill("NJets_{0}^{=3} & MET#geq 150",1);}
		    }
		  if(njets==4)
		    {
		      if(MET>=100 && MET<150)
			{ ll_2->Fill("NJets_{0}^{=4} & 100<MET<150",1);}
		      if(MET>=150)
			{ ll_2->Fill("NJets_{0}^{=4} & MET#geq 150",1);}
		    }
		}
	      if(goodjets.size()>=5 && goodjets.size()<=6)
		{
		  ll->Fill("NJets_{=0}^{5-6}",1);
		  if(MET>=100 && MET<150)
		    { ll_2->Fill("NJets_{0}^{5-6} & 100<MET<150",1);}
		  if(MET>=150)
		    { ll_2->Fill("NJets_{0}^{5-6} & MET#geq 150",1);}
		}
	      if(goodjets.size()>=7)
		{
		  ll->Fill("NJets_{=0}^{#geq 7}",1);
		  if(MET>=100 && MET<150)
		    { ll_2->Fill("NJets_{0}^{#geq7} & 100<MET<150",1);}
		  if(MET>=150)
		    { ll_2->Fill("NJets_{0}^{#geq7} & MET#geq 150",1);}
		}

	    }
	  if(BTags>=1)
	    { if(goodjets.size()>=2 && goodjets.size()<=4)
		{
		  ll->Fill("NJets_{#geq 1}^{2-4}",1);
		  if(MET>=100 && MET<150)
		    { ll_2->Fill("NJets_{#geq 1}^{2-4} & 100<MET<150",1);}
		  if(MET>=150)
		    { ll_2->Fill("NJets_{#geq 1}^{2-4} & MET#geq 150",1);}
		}
	      if(goodjets.size()>=5 && goodjets.size()<=6)
		{
		  ll->Fill("NJets_{#geq 1}^{5-6}",1);
		  if(MET>=100 && MET<150)
		    { ll_2->Fill("NJets_{#geq 1}^{5-6} & 100<MET<150",1);}
		  if(MET>=150)
		    { ll_2->Fill("NJets_{#geq 1}^{5-6} & MET#geq 150",1);}
		}
	      if(goodjets.size()>=7)
		{
		  ll->Fill("NJets_{#geq 1}^{#geq 7}",1);
		  if(MET>=100 && MET<150)
		    { ll_2->Fill("NJets_{#geq 1}^{#geq 7} & 100<MET<150",1);}
		  if(MET>=150)
		    { ll_2->Fill("NJets_{#geq 1}^{#geq 7} & MET#geq 150",1);}
		}
	      
	    }	
	  
	}			
      


      
    }// loop over entries

  
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
