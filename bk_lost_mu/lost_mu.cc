#define lost_mu_cxx
#include "lost_mu.h"
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

  lost_mu ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void lost_mu::EventLoop(const char *data,const char *inputFileList)
{ 
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  TString s_data=data;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;

  int survived_events =0,survived_vetohad=0,not_accepted=0,survived_accept=0,check=0,survived_failed_id=0,survived_failed_iso=0,survived_all=0,events_cr=0,events_id=0,events_failiso=0;
  int event_isol2=0,event_miniiso=0,event_failminiiso=0,event_failiso2=0,event_failiso=0,event_isol=0;
  int e_mu_event=0,events_gen_mu=0,events_hadtau=0;
  
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

      ///////////////////////   Event weight   ////////////////
      Double_t wt = 35.9*1000*(Weight);

      // filters
      if(strcmp(data,"2016")==1)
        {
          if(PrimaryVertexFilter!=1 || globalSuperTightHalo2016Filter!=1 || HBHEIsoNoiseFilter!=1 || HBHENoiseFilter!=1 || EcalDeadCellTriggerPrimiti\
veFilter!=1)    continue;
          if(BadPFMuonFilter == false) continue;
        }
      if(strcmp(data,"2017")==1)
        {
          if(PrimaryVertexFilter!=1 || globalSuperTightHalo2016Filter!=1 || HBHEIsoNoiseFilter!=1 || HBHENoiseFilter!=1 || EcalDeadCellTriggerPrimiti\
veFilter!=1)    continue;
          if(BadPFMuonFilter == false) continue;
        }

      
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

      // if(goodphoton.Pt()<100) continue; // remove events with 0 photons or pt<100
      // if(MET<100) continue;
      
      
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
      // if(dphi1>0.3 && dphi2>0.3) continue;
      // if((st<800 || goodphoton.Pt()<100) && (st<500 || goodphoton.Pt()<190)) continue;
      // if(isoElectronTracks!=0 || isoPionTracks!=0) continue;

      double mt_mu=0,mt_pho=0,mt_mupho=0;
      if(Muons->size()==1){
	mt_mu=sqrt(2*(*Muons)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*Muons)[0].Phi()))));
	if(mt_mu>100) continue;
	if( ((*Muons)[0].Pt() < 10) || abs((*Muons)[0].Eta()) > 2.5 ) continue;
      }
      mt_pho=sqrt(2*goodphoton.Pt()*MET*(1-cos(DeltaPhi(METPhi,goodphoton.Phi()))));
      if(Muons->size()==1) mt_mupho=sqrt(2*(goodphoton+(*Muons)[0]).Pt()*MET*(1-cos(DeltaPhi(METPhi,(goodphoton+(*Muons)[0]).Phi()))));   // MHT cut defined in the paper for isolated tracks

    
      

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

      nel[0]->Fill(NElectrons,wt);
      nmu[0]->Fill(NMuons,wt);


      ////////         Preselection cut                 //////////////////
      if((goodphoton.Pt()>100)&&(MET>=100)&&(dphi1>0.3)&&(dphi2>0.3)&&(((st>800 && goodphoton.Pt()>100) || (st>500 && goodphoton.Pt()>190)))&&(isoElectronTracks==0)&&(isoPionTracks==0)&&(goodjets.size()>=2)&&(NElectrons==0))
	{ 
	  //	  cout<<"Not here"<<endl;

	  survived_events+=1;

      
	  // Getting the Gen Information
	  
	  // ----------------------  Gen level -------------------------
	  vector<TLorentzVector> gen_el, gen_mu,gen_tau,gen_tau_el,gen_tau_mu,gen_ph;
	  
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
		  if((abs((*GenParticles_PdgId)[i])==15) && (abs((*GenParticles_ParentId)[i])==24))
		    { 
		      gen_tau.push_back((*GenParticles)[i]);
		    }
		  if((abs((*GenParticles_PdgId)[i])==22))
		    {
		      gen_ph.push_back((*GenParticles)[i]);
		    }		        
		}
	    }
	  	  
	  if(gen_el.size()==0 && gen_mu.size()==0 && gen_tau.size()==0) continue; // rejecting hadronic decays


	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////    Hadronic Tau events        ////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  bool hadtau = false;
	  if((NMuons==0)&&(NElectrons==0)&&(isoMuonTracks==0)&&(isoElectronTracks==0)&&(isoPionTracks==0))
	    {
	      if(gen_tau_el.size()==0 && gen_tau_mu.size()==0 && gen_tau.size()!=0)
		{ for(int i=0;i<gen_tau.size();i++)
		    { if(MinDr(gen_tau[i],goodjets)<0.3)
			{ events_hadtau+=1;
			  fill_hist(hadtau1,hadtau2,BTags,MET,goodjets.size(),wt);
			  hadtau=true;
			}
		    }
		}
	    }
	  
	  ///////////////////////////           back to lost mu case                  ///////////////////////////////////////////////////////////
	  
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

	  gen_el_size[0]->Fill(gen_el.size(),wt);
	  for(int i=0;i<gen_el.size();i++)
	    { gen_el_pt->Fill(gen_el[i].Pt(),wt);
	      gen_el_eta->Fill(gen_el[i].Eta(),wt);
	    }
      
	  gen_mu_size[0]->Fill(gen_mu.size(),wt);
	  for(int i=0;i<gen_mu.size();i++)
	    { gen_mu_pt->Fill(gen_mu[i].Pt(),wt);
	      gen_mu_eta->Fill(gen_mu[i].Eta(),wt);
	    }

	  gen_ph_size[0]->Fill(gen_ph.size(),wt);
	  for(int i=0;i<gen_ph.size();i++)
	    { gen_ph_pt->Fill(gen_ph[i].Pt(),wt);
	      gen_ph_eta->Fill(gen_ph[i].Eta(),wt);
	    }

	  survived_vetohad+=1;
	  nel[1]->Fill(NElectrons,wt);
	  nmu[1]->Fill(NMuons,wt);
	  // if(gen_mu.size()!=1) continue;
	  // if(gen_el.size()!=0) continue;
	  // if(NMuons>=2) continue;

	  if(gen_mu.size()!=0)
	    {
 
	      for(int i=0;i<gen_mu.size();i++)
		{
		  fill_hist(total1,total2,BTags,MET,goodjets.size(),wt);
      
      
		  // Out of acceptance ////////////////////////////////////////////////
      
		  bool acceptance = true;
		   if((abs(gen_mu[i].Eta())>2.5)||(gen_mu[i].Pt()<10))
		     { not_accepted+=1;
		       fill_hist(fail_accept1,fail_accept2,BTags,MET,goodjets.size(),wt);
		       acceptance =false;
		     }
		    
      
		  if(acceptance)
		    {
		      //////////////////////////////////////////////////////////////
		      survived_accept+=1;
		      nel[2]->Fill(NElectrons,wt);
		      nmu[2]->Fill(NMuons,wt);
		      gen_el_size[1]->Fill(gen_el.size(),wt);
		      gen_mu_size[1]->Fill(gen_mu.size(),wt);
		      gen_ph_size[1]->Fill(gen_ph.size(),wt);
     
      
      
		      // Failed Id  ////////////////////////////////////////////////

		      bool identified = true;
		      mindr_gen_rec_mu ->Fill(MinDr(gen_mu[i],*Muons),wt);
		      if((MinDr(gen_mu[i],*Muons)>0.2))
			{ events_id+=1;
			  fill_hist(fail_id1,fail_id2,BTags,MET,goodjets.size(),wt);
			  identified = false;
			}
			
		      if(identified)
			{
			  ////////////////////////////////////////////////////////////
			  survived_failed_id+=1;
			  nel[3]->Fill(NElectrons,wt);
			  nmu[3]->Fill(NMuons,wt);
			  gen_el_size[2]->Fill(gen_el.size(),wt);
			  gen_mu_size[2]->Fill(gen_mu.size(),wt);
			  gen_ph_size[2]->Fill(gen_ph.size(),wt);
      
      

			  //   Failed Iso ///////////////////////////////////////////
      
			  bool isol= true;
			  if(NMuons == 0)
			    { events_failiso+=1; 
			      fill_hist(fail_iso1,fail_iso2,BTags,MET,goodjets.size(),wt);
			      isol = false;
			    }

			  if(Muons->size()!=0 && !Muons_passIso)
			    { event_failiso2+=1;	}
			  if(NMuons>0)
			    { event_isol+=1;} 
			  if(Muons->size()!=0 && Muons_passIso)
			    { event_isol2+=1;	}
  

			  bool miniIso=false;
			  for(int i=0;i<Muons_MiniIso->size();i++)
			    { if((*Muons_MiniIso)[i]<0.1)
				{ miniIso = true;}
			    }
			  if(miniIso)
			    { event_miniiso+=1;}
			  else if(!miniIso)
			    { event_failminiiso+=1;}
      

      
			  if(isol)
			    {
      
      
			      /////////////////////////////////////////////////////////
			      survived_failed_iso+=1;
			      nel[4]->Fill(NElectrons,wt);
			      nmu[4]->Fill(NMuons,wt);
			      gen_el_size[3]->Fill(gen_el.size(),wt);
			      gen_mu_size[3]->Fill(gen_mu.size(),wt);
			      gen_ph_size[3]->Fill(gen_ph.size(),wt);
			      // 1 Lepton CR ///////////////////////////////////////////
			      bool cr = true;
			      if(NMuons == 1)
				{ events_cr+=1;
				  fill_hist(one_lep_cr1,one_lep_cr2,BTags,MET,goodjets.size(),wt);
				  cr = false;
				}
			      if(cr)
				{

				  /////////////////////////////////////////////////////////
				  survived_all+=1;
				  nel[5]->Fill(NElectrons,wt);
				  nmu[5]->Fill(NMuons,wt);
				  gen_el_size[4]->Fill(gen_el.size(),wt);
				  gen_mu_size[4]->Fill(gen_mu.size(),wt);
				  gen_ph_size[4]->Fill(gen_ph.size(),wt);
				} //cr
			    } // iso
			} // id
		    } // acceptance
		}//loop over all the gen muons
	
	    }//if(gen_mu.size()!=0)
	}//preselection 
    }// loop over entries
    
  cout<<"Events survived preselection         = "<<survived_events<<endl;
  cout<<"Events survived veto had             = "<<survived_vetohad<<endl;
  cout<<"Events hadronic tau                  = "<<events_hadtau<<endl;
  cout<<"Events not accepted by the detector  = "<<not_accepted<<endl;
  cout<<"Events survived acceptance cut       = "<<survived_accept<<endl;
  //cout<<" CHECK                               = "<<check<<endl;
  cout<<"Events failed ID                     = "<<events_id<<endl;
  cout<<"Events survived id cut               = "<<survived_failed_id<<endl;
  cout<<"Events failed Iso                    = "<<events_failiso<<endl;
  cout<<"Events failed pass iso (!passiso)    = "<<event_failiso2<<endl;
  cout<<"Events failed mini iso (MiniIso>0.2) ="<<event_failminiiso<<endl;
  cout<<"Events pass NElectrons>0             ="<<event_isol<<endl;
  cout<<"Events pass (passIso)                = "<<event_isol2<<endl;
  cout<<"Events pass miniiso <0.1             = "<<event_miniiso<<endl;
  cout<<"Events survived iso cut              = "<<survived_failed_iso<<endl;
  cout<<"Events in 1 lepton CR                = "<<events_cr<<endl;
  cout<<"Events survived all cut              = "<<survived_all<<endl;
  
  
}
	
  
bool lost_mu::electron_match_photon(TLorentzVector photon)
{ for(int i=0;i<Electrons->size();i++)
    { if(photon.DeltaR((*Electrons)[i])<0.2)
	{ return true;	}
      else
	{ return false;}
    }
}


void lost_mu::fill_hist(TH1D *h1, TH1D *h2,int btags,double met,int njets,double wt)
{
  //h1 = TH1D("h1","lost electron in b-jets and njets bins",6,1,7);
  //h2 = TH1D("h2","lost electron in all the bins",16,1,17);
  
  h1->Fill("NJets_{=0}^{2-4}",0);
  h1->Fill("NJets_{#geq 1}^{2-4}",0);
  //  h1->Fill("NJets_{#geq 2}^{2-4}",0);
  h1->Fill("NJets_{=0}^{5-6}",0);
  h1->Fill("NJets_{#geq 1}^{5-6}",0);
  // h1->Fill("NJets_{#geq 2}^{5-6}",0);
  h1->Fill("NJets_{=0}^{#geq 7}",0);
  h1->Fill("NJets_{#geq 1}^{#geq 7}",0);
  // h1->Fill("NJets_{#geq 2}^{#geq 7}",0);
  
  //int njets = goodjets.size();
  
  h2->Fill("NJets_{0}^{=2} & 100<MET<150",0);
  h2->Fill("NJets_{0}^{=2} & MET#geq 150",0);
  h2->Fill("NJets_{0}^{=3} & 100<MET<150",0);
  h2->Fill("NJets_{0}^{=3} & MET#geq 150",0);
  h2->Fill("NJets_{0}^{=4} & 100<MET<150",0);
  h2->Fill("NJets_{0}^{=4} & MET#geq 150",0);
  h2->Fill("NJets_{0}^{5-6} & 100<MET<150",0);
  h2->Fill("NJets_{0}^{5-6} & MET#geq 150",0);
  h2->Fill("NJets_{0}^{#geq7} & 100<MET<150",0);
  h2->Fill("NJets_{0}^{#geq7} & MET#geq 150",0);
  h2->Fill("NJets_{#geq 1}^{2-4} & 100<MET<150",0);
  h2->Fill("NJets_{#geq 1}^{2-4} & MET#geq 150",0);
  h2->Fill("NJets_{#geq 1}^{5-6} & 100<MET<150",0);
  h2->Fill("NJets_{#geq 1}^{5-6} & MET#geq 150",0);
  h2->Fill("NJets_{#geq 1}^{#geq 7} & 100<MET<150",0);
  h2->Fill("NJets_{#geq 1}^{#geq 7} & MET#geq 150",0);
  // h2->Fill("NJets_{#geq 2}^{2-4} & 100<MET<150",0);
  // h2->Fill("NJets_{#geq 2}^{2-4} & MET#geq 150",0);
  // h2->Fill("NJets_{#geq 2}^{5-6} & 100<MET<150",0);
  // h2->Fill("NJets_{#geq 2}^{5-6} & MET#geq 150",0);
  // h2->Fill("NJets_{#geq 2}^{#geq 7} & 100<MET<150",0);
  // h2->Fill("NJets_{#geq 2}^{#geq 7} & MET#geq 150",0);
  
  if(btags==0)
    {
      if(njets>=2 && njets<=4)
	{		  
	  h1->Fill("NJets_{=0}^{2-4}",wt);
	  if(njets==2)
	    {  
	      if(met>=100 && met<150)
		{ h2->Fill("NJets_{0}^{=2} & 100<MET<150",wt);}
	      if(met>=150)
		{ h2->Fill("NJets_{0}^{=2} & MET#geq 150",wt);}
	    }
	  if(njets==3)
	    {
	      if(met>=100 && met<150)
		{ h2->Fill("NJets_{0}^{=3} & 100<MET<150",wt);}
	      if(met>=150)
		{ h2->Fill("NJets_{0}^{=3} & MET#geq 150",wt);}
	    }
	  if(njets==4)
	    {
	      if(met>=100 && met<150)
		{ h2->Fill("NJets_{0}^{=4} & 100<MET<150",wt);}
	      if(met>=150)
		{ h2->Fill("NJets_{0}^{=4} & MET#geq 150",wt);}
	    }
	}
      if(njets>=5 && njets<=6)
	{
	  h1->Fill("NJets_{=0}^{5-6}",wt);
	  if(met>=100 && met<150)
	    { h2->Fill("NJets_{0}^{5-6} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2->Fill("NJets_{0}^{5-6} & MET#geq 150",wt);}
	}
      if(njets>=7)
	{
	  h1->Fill("NJets_{=0}^{#geq 7}",wt);
	  if(met>=100 && met<150)
	    { h2->Fill("NJets_{0}^{#geq7} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2->Fill("NJets_{0}^{#geq7} & MET#geq 150",wt);}
	}

    }
  if(btags>=1)
    { if(njets>=2 && njets<=4)
	{
	  h1->Fill("NJets_{#geq 1}^{2-4}",wt);
	  if(met>=100 && met<150)
	    { h2->Fill("NJets_{#geq 1}^{2-4} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2->Fill("NJets_{#geq 1}^{2-4} & MET#geq 150",wt);}
	}
      if(njets>=5 && njets<=6)
	{
	  h1->Fill("NJets_{#geq 1}^{5-6}",wt);
	  if(met>=100 && met<150)
	    { h2->Fill("NJets_{#geq 1}^{5-6} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2->Fill("NJets_{#geq 1}^{5-6} & MET#geq 150",wt);}
	}
      if(njets>=7)
	{
	  h1->Fill("NJets_{#geq 1}^{#geq 7}",wt);
	  if(met>=100 && met<150)
	    { h2->Fill("NJets_{#geq 1}^{#geq 7} & 100<MET<150",wt);}
	  if(met>=150)
	    { h2->Fill("NJets_{#geq 1}^{#geq 7} & MET#geq 150",wt);}
	}
	      
    }
  // if(btags>=2)
  //   { if(njets>=2 && njets<=4)
  // 	{
  // 	  h1->Fill("NJets_{#geq 2}^{2-4}",wt);
  // 	  if(met>=100 && met<150)
  // 	    { h2->Fill("NJets_{#geq 2}^{2-4} & 100<MET<150",wt);}
  // 	  if(met>=150)
  // 	    { h2->Fill("NJets_{#geq 2}^{2-4} & MET#geq 150",wt);}
  // 	}
  //     if(njets>=5 && njets<=6)
  // 	{
  // 	  h1->Fill("NJets_{#geq 2}^{5-6}",wt);
  // 	  if(met>=100 && met<150)
  // 	    { h2->Fill("NJets_{#geq 2}^{5-6} & 100<MET<150",wt);}
  // 	  if(met>=150)
  // 	    { h2->Fill("NJets_{#geq 2}^{5-6} & MET#geq 150",wt);}
  // 	}
  //     if(njets>=7)
  // 	{
  // 	  h1->Fill("NJets_{#geq 2}^{#geq 7}",wt);
  // 	  if(met>=100 && met<150)
  // 	    { h2->Fill("NJets_{#geq 2}^{#geq 7} & 100<MET<150",wt);}
  // 	  if(met>=150)
  // 	    { h2->Fill("NJets_{#geq 2}^{#geq 7} & MET#geq 150",wt);}
  // 	}
	      
  //   }
}





