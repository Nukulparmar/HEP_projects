#define fake_photon_cxx
#include "fake_photon.h"
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

  fake_photon ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void fake_photon::EventLoop(const char *data,const char *inputFileList)
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

      if(Electrons->size()>1 || Muons->size()!=0) continue;

      bool is_electron=false,is_photon=false;
      if(Electrons->size()==1)
	{ if((*Electrons)[0].Pt()>100)
	    is_electron = true;             // we think that electron is faking as photon, so we select events with one electron which has photon like cuts.
	}

      double pt = 0.0;
      int goodphoton_index = -1;
      TLorentzVector goodphoton;
      for(int i=0;i<Photons->size();i++)
	{ if(((*Photons)[i].Pt()>100)&&(TMath::Abs((*Photons)[i].Eta())<2.4)&&((*Photons_fullID)[i]==1))
	    { goodphotons.push_back((*Photons)[i]);
	      if(pt<(*Photons)[i].Pt())
		{ pt = (*Photons)[i].Pt();
		  goodphoton_index = i;
		}
	      
	    }
	  
	  
	}
      sortTLorVec(&goodphotons);
      if(goodphotons.size()!=0)
	{goodphoton = goodphotons[0];}

      
      bool matched = electron_match_photon(goodphoton); // checking reco photon and reco electron vectors
      bool goodphoton_hasPixelSeed=true;
      if(goodphoton_index>=0)
	{ if((*Photons_hasPixelSeed)[goodphoton_index]<0.001) goodphoton_hasPixelSeed=false;
	  if(!goodphoton_hasPixelSeed && goodphoton.Pt()>100 && !matched) is_photon = true;
	}

      TLorentzVector object;
      bool object_iselectron;
      if(is_electron && !is_photon)
	{ object = (*Electrons)[0];
	  object_iselectron = true;
	}
      else if(!is_electron && is_photon)
	{ object = (*Photons)[0];
	  object_iselectron = false;
	}
      else continue;
      //checking if gen electron matches the reco photon
      bool fakephoton = false;
      if(!object_iselectron)                          
	{ double dr_ph_genobj=100;
	  for(int i=0;i<GenParticles->size();i++)
	    { if((*GenParticles)[i].Pt()!=0)
		{ double dr = goodphoton.DeltaR((*GenParticles)[i]);
		  if(dr<0.2 && (abs((*GenParticles_PdgId)[i])==11) && (abs((*GenParticles_ParentId)[i])<=24))
		    { fakephoton = true;}
		
		  if(dr<dr_ph_genobj)
		    { dr_ph_genobj = dr;}
		}
	    }

	  for(int i=0;i<GenParticles->size();i++)
	    {
	      if((*GenParticles)[i].Pt()!=0)
		{ double dr = goodphoton.DeltaR((*GenParticles)[i]);
		  if(dr<0.2 && (abs((*GenParticles_PdgId)[i])==22) && (((*GenParticles)[i].Pt()/goodphoton.Pt())>0.95) && (((*GenParticles)[i].Pt()/goodphoton.Pt())>1.05))
		    { fakephoton = false;}
		}
	    }
	  if(!fakephoton) continue;
	  if(isoMuonTracks!=0 || isoPionTracks!=0 || isoElectronTracks!=0) continue;
	}
      

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

      

      if(MET>100 && goodjets.size()>=2 && (dphi1>0.3 && dphi2 >0.3) && ht>100 && goodphoton.Pt()>100)
	{ survived_events+=1;

	  if(BTags==0)
	    { if(goodjets.size()>=2 && goodjets.size()<=4)
		{		  
		  faking_photon->Fill(1,1);
		}
	      if(goodjets.size()>=5 && goodjets.size()<=6)
		{
		  faking_photon->Fill(2,1);
		}
	      if(goodjets.size()>=7)
		{
		  faking_photon->Fill(3,1);
		}

	    }
	  if(BTags>=1)
	    { if(goodjets.size()>=2 && goodjets.size()<=4)
		{
		  faking_photon->Fill(4,1);
		}
	      if(goodjets.size()>=5 && goodjets.size()<=6)
		{
		  faking_photon->Fill(5,1);
		}
	      if(goodjets.size()>=7)
		{
		  faking_photon->Fill(6,1);
		}
	      
	    }	
	  
	}			
      


      
    }// loop over entries

  
}
	
  
bool fake_photon::electron_match_photon(TLorentzVector photon)
{ for(int i=0;i<Electrons->size();i++)
    { if(photon.DeltaR((*Electrons)[i])<0.2)
	{ return true;
	}
    
      else
	{ return false;}

    }
}
