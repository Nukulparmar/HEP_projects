#define had_tau_cxx
#include "had_tau.h"
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

  had_tau ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void had_tau::EventLoop(const char *data,const char *inputFileList)
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
      if(electron_match_photon(goodphoton)) continue; // don't want electron matched to good photon

      if(Muons->size()>0) continue; // veto muons
      if(Electrons->size()>0) continue; //  no reco electrons
      if(isoMuonTracks!=0 || isoPionTracks!=0 || isoElectronTracks!=0) continue;

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
      int nGenMu=0,nGenEle=0,nGenTau=0,nGenMuFmTau=0,nGenEleFmTau=0;

      for(int i=0;i<GenParticles->size();i++)
	{ if((*GenParticles)[i].Pt()!=0)
	    {
	      if((abs((*GenParticles_PdgId)[i])==11) && (abs((*GenParticles_ParentId)[i])<=24) && ((*GenParticles_Status)[i]==1))
		{ 
		  gen_el.push_back((*GenParticles)[i]);
		  if(abs((*GenParticles_ParentId)[i])==15)
		    {gen_tau_el.push_back((*GenParticles)[i]);}
		}
	      if((abs((*GenParticles_PdgId)[i])==13) && (abs((*GenParticles_ParentId)[i])<=24) && ((*GenParticles_Status)[i]==1))
		{ 
		  gen_mu.push_back((*GenParticles)[i]);
		  if(abs((*GenParticles_ParentId)[i])==15)
		    {gen_tau_mu.push_back((*GenParticles)[i]);}
		}
	      if((abs((*GenParticles_PdgId)[i])==15) && (abs((*GenParticles_ParentId)[i])<=24) && ((*GenParticles_Status)[i]==1))
		{ 
		  gen_tau.push_back((*GenParticles)[i]);
		}
	      
	    }
	}

      if(gen_el.size()!=0 && gen_mu.size()!=0) continue; // only hadronic decays
      if(gen_tau_el.size()!=0 && gen_tau_mu.size()!=0) continue;
      
      sortTLorVec(&gen_tau);

      // checking if the photon is real or fake

      bool realphoton = true;
      int match_el=0,match_p=0;
      double mindr_ph_genobj = 100,match_el_pt=0.0,match_ph_pt=0.0;
      for(int i=0;i<GenParticles->size();i++)
	{ if((*GenParticles)[i].Pt()!=0)
	    { double dr1=goodphoton.DeltaR((*GenParticles)[i]);
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

      if(match_el==1 && match_p==0) realphoton = false;
      
      if(!realphoton) continue;
      

      // tau matching with jets 

      double mindr_tau_jet = 100.0;
      int index_jet_match_tau = -1;
      bool tau_match_jet = false;
      
	 for(int j=0;j<gen_tau.size();j++)
	    {
	      { for(int i=0;i<Jets->size();i++)
		  {
		    if((*Jets)[i].Pt()>30 && abs((*Jets)[i].Eta())<=2.5)
		      {
			double dr = gen_tau[j].DeltaR((*Jets)[i]);
			if(dr<mindr_tau_jet)
			  { mindr_tau_jet = dr;
			    index_jet_match_el = i;
			  }
		      }
		  }
	      }
	    

	  // if(index_jet_match_el>=0 && ((*Jets)[index_jet_match_el].Pt())/((*Electrons)[0].Pt())<1.0) continue;
	  // if(index_jet_match_el<0) continue;
	}
      // if(gamma_matching_jet_index>=0 && ((*Jets)[gamma_matching_jet_index].Pt())/((*Electrons)[0].Pt())<1.0) continue;
      // if(gamma_matching_jet_index<0) continue;


      if(MET>100 && goodjets.size()>=2 && (dphi1>0.3 && dphi2 >0.3) && ht>100 && goodphoton.Pt()>100)
	{ survived_events+=1;

	  if(BTags==0)
	    { if(goodjets.size()>=2 && goodjets.size()<=4)
		{		  
		  hadtau->Fill(1,1);
		}
	      if(goodjets.size()>=5 && goodjets.size()<=6)
		{
		  hadtau->Fill(2,1);
		}
	      if(goodjets.size()>=7)
		{
		  hadtau->Fill(3,1);
		}

	    }
	  if(BTags>=1)
	    { if(goodjets.size()>=2 && goodjets.size()<=4)
		{
		  hadtau->Fill(4,1);
		}
	      if(goodjets.size()>=5 && goodjets.size()<=6)
		{
		  hadtau->Fill(5,1);
		}
	      if(goodjets.size()>=7)
		{
		  hadtau->Fill(6,1);
		}
	      
	    }	
	  
	}			
      


      
    }// loop over entries
  cout<<"Events survived = "<<survived_events;
  
}
	
  
bool had_tau::electron_match_photon(TLorentzVector photon)
{ for(int i=0;i<Electrons->size();i++)
    { if(photon.DeltaR((*Electrons)[i])<0.2)
	{ return true;
	}
    
      else
	{ return false;}

    }
}
