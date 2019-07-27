#define AnalyzeLightBSM_cxx
#include "AnalyzeLightBSM.h"
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

  AnalyzeLightBSM ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList)
{ 
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  TString s_data=data;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;

  int count=0;
  //////////////////////////////////////
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

      ///////////////         GEN LEVEL       ////////////////////////

      
      vector<TLorentzVector> gen_photon,gen_z,gen_w,nu,gen_w_el,gen_w_mu,gen_w_lep,gen_w_tau;

      for(int i=0;i<GenParticles->size();i++)
	{
	  if((*GenParticles_PdgId)[i]==22)
	    { gen_photon.push_back((*GenParticles)[i]);
	      h_gen_photon_pt->Fill((*GenParticles)[i].Pt(),wt);
	    }

	  if((*GenParticles_PdgId)[i]==23)
	    {	      gen_z.push_back((*GenParticles)[i]);	    }

	  if(abs((*GenParticles_PdgId)[i])==24)
	    {	      gen_w.push_back((*GenParticles)[i]);	    }

	  if((abs((*GenParticles_PdgId)[i])==12)||(abs((*GenParticles_PdgId)[i])==14)||(abs((*GenParticles_PdgId)[i])==16))
	    { nu.push_back((*GenParticles)[i]);}

	  if((abs((*GenParticles_PdgId)[i])==11)&&(abs((*GenParticles_ParentId)[i])==24))
	    { gen_w_el.push_back((*GenParticles)[i]);
	      gen_w_lep.push_back((*GenParticles)[i]);
	    }

	  if((abs((*GenParticles_PdgId)[i])==13)&&(abs((*GenParticles_ParentId)[i])==24))
	    { gen_w_mu.push_back((*GenParticles)[i]);
	      gen_w_lep.push_back((*GenParticles)[i]);
	    }

	  if((abs((*GenParticles_PdgId)[i])==15)&&(abs((*GenParticles_ParentId)[i])==24))
	    { gen_w_tau.push_back((*GenParticles)[i]);
	    }
  
	}

      
      sortTLorVec(&gen_photon);

      if(gen_photon.size()!=0)
	{ h_gen_lead_photon_pt->Fill(gen_photon[0].Pt(),wt);}


      ////////////////////////////  For QCD ->  Finding the miss-measurement contributing to MET //////////////////////

      double dr = 1000, eff =0.0; int index = -1;
      for(int i=0;i<GenJets->size();i++)
	{
	  dr = 1000.0;
	  index = -1;
	  for(int j=0;j<Jets->size();j++)
	    {
	      
	      if(dr>(*GenJets)[i].DeltaR((*Jets)[j]))
		{ dr = (*GenJets)[i].DeltaR((*Jets)[j]);
		  index = j;
		}
	    }
	  if(dr<0.2)
	    {
	      eff = fabs(((*GenJets)[i].Pt()-(*Jets)[index].Pt())/(*GenJets)[i].Pt());
	      h_eff->Fill(eff,wt);
	      eff_vs_jetpt->Fill((*GenJets)[i].Pt(),eff,wt);
	    }
	}
      
      ////////////////      RECO LEVEL     //////////////////////////////


      if(Photons->size()!=0)
	{ h_lead_photon_pt->Fill((*Photons)[0].Pt(),wt);}

      for(int i=0;i<Photons->size();i++)
	{ h_photon_pt->Fill((*Photons)[i].Pt(),wt);
	  h_photon_eta->Fill((*Photons)[i].Eta(),wt);
	  h_photon_phi->Fill((*Photons)[i].Phi(),wt);
	}

      double ht = 0;
      vector<TLorentzVector> goodjets,btag_jets;
      for(int i=0;i<Jets->size();i++)
	{ if(((*Jets)[i].Pt()>30)&&(TMath::Abs((*Jets)[i].Eta())<3))
	    { ht += (*Jets)[i].Pt();
	      goodjets.push_back((*Jets)[i]);
	      
	    }
	}

      h_HT->Fill(ht,wt);
      h_Njets->Fill(goodjets.size(),wt);
      h_Nbjets->Fill(BTags,wt);
      sortTLorVec(&goodjets);

      h_MET->Fill(MET,wt);
      h_MHT->Fill(MHT,wt);

      h_isoeltracks->Fill(isoElectronTracks,wt);
      h_isomutracks->Fill(isoMuonTracks,wt);
      h_isopitracks->Fill(isoPionTracks,wt);


      ///  min dr between photon and jet
       dr = 100.0;
       index = -1;
      for(int i=0;i<Photons->size();i++)
	{
	  dr = 100.0;
	  index =-1;
	  for(int j=0;j<goodjets.size();j++)
	    { if(dr>=goodjets[j].DeltaR((*Photons)[i]))
		{ dr = goodjets[j].DeltaR((*Photons)[i]);
		  index = j;
		}
	      if(index!=-1)
		{ h_dr_photon_jets->Fill(dr,wt);}
	    }
	}

      for(int i=0;i<Electrons->size();i++)
	{ h_el_pt->Fill((*Electrons)[i].Pt(),wt);}
      for(int i=0;i<Muons->size();i++)
	{ h_mu_pt->Fill((*Muons)[i].Pt(),wt);}

      

      //  Electron Faking Photons

      if((GenElectrons->size()>0)&&(Photons->size()>0))
	{ for(int i=0;i<GenElectrons->size();i++)
	    { dr = 1000.0;
	      index = -1;
	      for(int j=0;j<Photons->size();j++)
		{ if(dr>(*GenElectrons)[i].DeltaR((*Photons)[j]))
		    { dr = (*GenElectrons)[i].DeltaR((*Photons)[j]);
		      index = j;
		    }
		}
	      if(index!=-1)
		{ h_dr_el_fake_photon->Fill(dr,wt);}
	      if(dr<0.1)
		{ genel_vs_recphoton->Fill((*GenElectrons)[i].Pt(),(*Photons)[index].Pt());
		  el_ph_fake_corr->Fill((*GenElectrons)[i].Pt(),(*Photons)[index].Pt()/(*GenElectrons)[i].Pt(),wt);
		}
	      
	    }
	}
      
      //   Jet faking electron

       for(int i=0;i<GenJets->size();i++)
	    { dr = 1000.0;
	      index = -1;
	      for(int j=0;j<Electrons->size();j++)
		{ if(dr>(*GenJets)[i].DeltaR((*Electrons)[j]))
		    { dr = (*GenJets)[i].DeltaR((*Electrons)[j]);
		      index = j;
		    }
		}
	      if(index!=-1)
		{ h_dr_jet_fake_el->Fill(dr,wt);}
	      if(dr<0.1)
		{ genjet_vs_recel->Fill((*GenJets)[i].Pt(),(*Electrons)[index].Pt());
		  jet_el_fake_corr->Fill((*GenJets)[i].Pt(),(*Electrons)[index].Pt()/(*GenJets)[i].Pt(),wt);
		}
	      
	    }

       //  Jet faking photons

       for(int i=0;i<GenJets->size();i++)
	    { dr = 1000.0;
	      index = -1;
	      for(int j=0;j<Photons->size();j++)
		{ if(dr>(*GenJets)[i].DeltaR((*Photons)[j]))
		    { dr = (*GenJets)[i].DeltaR((*Photons)[j]);
		      index = j;
		    }
		}
	      if(index!=-1)
		{ h_dr_jet_fake_ph->Fill(dr,wt);}
	      if(dr<0.1)
		{ genjet_vs_recph->Fill((*GenJets)[i].Pt(),(*Photons)[index].Pt());
		  jet_ph_fake_corr->Fill((*GenJets)[i].Pt(),(*Photons)[index].Pt()/(*GenJets)[i].Pt(),wt);
		}
	      
	    }
       


       
	  
      

	  
    }// loop over entries

  //cout<<"the analysed event number is = "<<count<<endl;
}
	
  
