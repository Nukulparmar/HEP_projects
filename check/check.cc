#define check_cxx
#include "check.h"
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

  check ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void check::EventLoop(const char *data,const char *inputFileList)
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
      
      num_events->Fill(wt);
      for(int i=0;i<Photons->size();i++)
	{
	  all_ph_pt->Fill((*Photons)[i].Pt(),wt);
	  all_ph_eta->Fill((*Photons)[i].Eta(),wt);
	  if(((*Photons)[i].Pt()>100)&&(TMath::Abs((*Photons)[i].Eta())<2.4)&&((*Photons_fullID)[i]==1)&&((*Photons_hasPixelSeed)[i]<0.000001))
	    { goodphotons.push_back((*Photons)[i]);
	      
	    }
	  
	}
      TLorentzVector goodphoton;
      sortTLorVec(&goodphotons);
      if(goodphotons.size()!=0)
	{goodphoton = goodphotons[0];
	  lead_ph_pt->Fill(goodphoton.Pt(),wt);
	  lead_ph_eta->Fill(abs(goodphoton.Eta()),wt);
	}
      bool temp=Photons_electronFakes;
      cout<<temp<<endl;
      ph_elfakes->Fill(boolvalue(Photons_electronFakes),wt);
      ph_fullId->Fill(boolvalue(Photons_fullID),wt);
      ph_pixelseed->Fill(boolvalue(Photons_hasPixelSeed),wt);
      for(int i=0;i<Photons_isEB->size();i++)
	{      ph_isEB->Fill((*Photons_isEB)[i],wt);     }
      for(int i=0;i<Photons_passElectronVeto->size();i++)
	{      ph_passelveto->Fill((*Photons_passElectronVeto)[i],wt);  }
      for(int i=0;i<Photons_pfGammaIso->size();i++)
	{      ph_pfgammaiso->Fill((*Photons_pfGammaIso)[i],wt);}            

      for(int i=0;i<Electrons->size();i++)
      	{ all_el_pt->Fill((*Electrons)[i].Pt(),wt);
      	  all_el_eta->Fill(abs((*Electrons)[i].Eta()),wt);
      	}

      for(int i=0;i<Muons->size();i++)
      	{ all_mu_pt->Fill((*Muons)[i].Pt(),wt);
      	  all_mu_eta->Fill(abs((*Muons)[i].Eta()),wt);
      	}

      // for(int i=0;i<Taus->size();i++)
      // 	{ all_mu_pt->Fill((*Taus)[i].Pt(),wt);
      // 	  all_mu_eta->Fill(abs((*Taus)[i].Eta()),wt);
      // 	}

      btags->Fill(BTags,wt);
      el_size->Fill(Electrons->size(),wt);
      el_mediumId->Fill(boolvalue(Electrons_mediumID),wt);
      for(int i=0;i<Electrons_MiniIso->size();i++)
	{el_miniIso->Fill((*Electrons_MiniIso)[i],wt);}
      el_passIso->Fill(boolvalue(Electrons_passIso),wt);
      el_tightId->Fill(boolvalue(Electrons_tightID),wt);
      el_isotrack->Fill(isoElectronTracks,wt);
      el_nelectrons->Fill(NElectrons,wt);



      mu_size->Fill(Muons->size(),wt);
      mu_mediumId->Fill(boolvalue(Muons_mediumID),wt);
      for(int i=0;i<Muons_MiniIso->size();i++)
	{      mu_miniIso->Fill((*Muons_MiniIso)[i],wt);}
      mu_passIso->Fill(boolvalue(Muons_passIso),wt);
      mu_tightId->Fill(boolvalue(Muons_tightID),wt);
      mu_isotrack->Fill(isoMuonTracks,wt);
      mu_nmuons->Fill(NMuons,wt);



      pion_isotrack->Fill(isoPionTracks,wt);
      
      

      
      

      
      // ----------------------  Gen level -------------------------
      vector<TLorentzVector> gen_el, gen_mu,gen_tau,gen_tau_el,gen_tau_mu,gen_ph;
     
      for(int i=0;i<GenParticles->size();i++)
      	{ if((*GenParticles)[i].Pt()!=0)
      	    {
      	      if((abs((*GenParticles_PdgId)[i])==11) && (abs((*GenParticles_ParentId)[i])==24) && ((*GenParticles_Status)[i]==1))
      	      	{ 
      	      	  gen_el.push_back((*GenParticles)[i]);
      	      	  if(abs((*GenParticles_ParentId)[i])==15)
      	      	    {gen_tau_el.push_back((*GenParticles)[i]);}
      	      	}
      	      if((abs((*GenParticles_PdgId)[i])==13) && (abs((*GenParticles_ParentId)[i])==24) && ((*GenParticles_Status)[i]==1))
      	      	{ 
      	      	  gen_mu.push_back((*GenParticles)[i]);
      	      	  if(abs((*GenParticles_ParentId)[i])==15)
      	      	    {gen_tau_mu.push_back((*GenParticles)[i]);}
      	      	}
      	      if((abs((*GenParticles_PdgId)[i])==15) && (abs((*GenParticles_ParentId)[i])==24) && ((*GenParticles_Status)[i]==1))
      	      	{ 
      	      	  gen_tau.push_back((*GenParticles)[i]);
      	      	}
      	      if((abs((*GenParticles_PdgId)[i]==22) && ((*GenParticles_Status)[i]==1)))
      	      	{ gen_ph.push_back((*GenParticles)[i]);}
	      
      	    }
      	  //cout<<i<<endl;
      	}
      
      if(gen_el.size()!=0)
	{ for(int i=0;i<gen_el.size();i++)
	    { gen_el_pt->Fill(gen_el[i].Pt(),wt);
	      gen_el_eta->Fill(abs(gen_el[i].Eta()),wt);
	    }
	}
      
      genel_size1->Fill(gen_el.size(),wt);
      genel_size2->Fill(GenElectrons->size(),wt);
      genel_fromtau->Fill(boolvalue(GenElectrons_fromTau),wt);
      if(gen_mu.size()!=0)
      	{
      	  for(int i=0;i<gen_mu.size();i++)
      	    { gen_mu_pt->Fill(gen_mu[i].Pt(),wt);
      	      gen_mu_eta->Fill(abs(gen_mu[i].Eta()),wt);
      	    }
      	}
      
      genmu_size1->Fill(gen_mu.size(),wt);
      genmu_size2->Fill(GenMuons->size(),wt);
      genmu_fromtau->Fill(boolvalue(GenMuons_fromTau),wt);
      
      for(int i=0;i<gen_tau.size();i++)
      	{ gen_tau_pt->Fill(gen_tau[i].Pt(),wt);
      	  gen_tau_eta->Fill(abs(gen_tau[i].Eta()),wt);
      	}
      gentau_size1->Fill(gen_tau.size(),wt);
      gentau_size2->Fill(GenTaus->size(),wt);
      gentau_had->Fill(boolvalue(GenTaus_had),wt);

      for(int i=0;i<GenTops->size();i++)
      	{ gentop_pt->Fill((*GenTops)[i].Pt(),wt);
      	  gentop_eta->Fill(abs((*GenTops)[i].Eta()),wt);
      	}
      gentop_size->Fill(GenTops->size(),wt);
      

     

      // ///   Additional plots

      double dr1 = 100,dr2=100;
      for(int i=0;i<Electrons->size();i++)
      	{ dr1 = mindr((*Electrons)[i],*Photons);
      	  mindr_recoel_recoph->Fill(dr1,wt);	  
      	  dr2 = mindr((*Electrons)[i],gen_el);
      	  mindr_recoel_genel->Fill(dr2,wt);
      	}
      dr1=100;dr2=100;
      for(int i=0;i<Photons->size();i++)
      	{ dr1 = mindr((*Photons)[i],gen_ph);
      	  mindr_recoph_genph->Fill(dr1,wt);
      	  dr2 = mindr((*Photons)[i],gen_el);
      	  mindr_recoph_genel->Fill(dr2,wt);
      	}
      //cout<<"Not Here"<<endl;
      
      for(int i=0;i<Electrons->size();i++)
      	{ for(int j=0;j<gen_el.size();j++)
      	    { genel_recoel_eta->Fill(abs(gen_el[j].Eta()),abs((*Electrons)[i].Eta()),wt);}
      	}
      for(int i=0;i<Photons->size();i++)
      	{ for(int j=0;j<gen_el.size();j++)
      	    { genel_recoph_eta->Fill(abs(gen_el[j].Eta()),abs((*Photons)[i].Eta()),wt);}
      	  for(int j=0;j<gen_ph.size();j++)
      	    { genph_recoph_eta->Fill(abs(gen_ph[j].Eta()),abs((*Photons)[i].Eta()),wt);}
      	  for(int j=0;j<Electrons->size();j++)
      	    { recoel_recoph_eta->Fill(abs((*Electrons)[j].Eta()),abs((*Photons)[i].Eta()),wt);}
      	}
	      
     
      
      
      
    }// loop over entries

}
	

bool check::electron_match_photon(TLorentzVector photon)
{ for(int i=0;i<Electrons->size();i++)
    { if(photon.DeltaR((*Electrons)[i])<0.2)
	{ return true;
	}
    
      else
	{ return false;}

    }
}

/*void check::fill_histo(TH1D hist1,TH1D hist2,vector<TLorentzVector> particle,double wt)
{ for(int i=0;i<particle.size();i++)
    { hist1.Fill(particle[i].Pt(),wt);
      hist2.Fill(abs(particle[i].Eta()),wt);
    }
    }*/

double check::mindr(TLorentzVector v1,vector<TLorentzVector> v2)
{
  double dr = 60;
  for(int j=0;j<v2.size();j++)
    { if(dr>=v1.DeltaR(v2[j]))
	{ dr = v1.DeltaR(v2[j]);}
    }
  return dr;
}

int check::boolvalue(bool a)
{
  if(a)
    { return 1;}
  if(!a)
    { return 0;}
}
      
