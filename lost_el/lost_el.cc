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
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" " << "year" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *year          = argv[4]; 
  lost_el ana(inputFileList, outFileName, data, year);
  cout << "dataset " << data << " " << " year " << year << endl;
  
  ana.EventLoop(data,inputFileList,year);

  return 0;
}

void lost_el::EventLoop(const char *data,const char *inputFileList, const char *year)
{ 
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;
  
  TString s_data=data;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
 
  //  pileup reweighting
  
  TFile *pufile;
  if(strcmp("2016",year)==0)
    {
      pufile = TFile::Open("/home/work/nparmar/SUSY_Photon_by_bhumi/SignalRegion/PileupHistograms_0121_69p2mb_pm4p6.root","READ");
    }
  if(strcmp(year,"2017")==0)
    {
      pufile = TFile::Open("/home/work/nparmar/pileup_files/PileupHistograms_0328_63mb_pm5.root","READ");
    }
  if(strcmp(year,"2018")==0)
    {
      pufile = TFile::Open("/home/work/nparmar/pileup_files/PileupHistograms_0118_63mb_pm5.root","READ");
    }
  //choose central, up, or down
  TH1* puhist = (TH1*)pufile->Get("pu_weights_down");

  
//   // Electron scalefactors
//   bool applyEGMSFs = 0;
//   TFile *f_EGMSF1=TFile::Open("scaleFactors.root");
//   TFile *f_EGMSF2=TFile::Open("egammaEffi.txt_EGM2D.root");
//   TH2F *h2_EGMSFV=(TH2F*)f_EGMSF1->Get("GsfElectronToCutBasedSpring15V");
//   TH2F *h2_EGMSFMVA=(TH2F*)f_EGMSF1->Get("GsfElectronToMVAVLooseTightIP2D");
//   TH2F *h2_EGMSFMiniIso=(TH2F*)f_EGMSF1->Get("MVAVLooseElectronToMini");
//   TH2F *h2_EGMSF1=(TH2F*)h2_EGMSFV->Clone("EGMSF1");
//   h2_EGMSF1->Multiply(h2_EGMSFMiniIso);// h2_EGMSF1->Divide(h2_EGMSFMVA);                                                                             
//   TH2F *h2_EGMSF2=(TH2F*)f_EGMSF2->Get("EGamma_SF2D");
//   cout<<"applying EGM SFs to electrons? "<<applyEGMSFs<<endl;
// bool applyEGMSFs = 0;
//   TFile *f_EGMSF1=TFile::Open("scaleFactors.root");
//   TFile *f_EGMSF2=TFile::Open("egammaEffi.txt_EGM2D.root");
//   TH2F *h2_EGMSFV=(TH2F*)f_EGMSF1->Get("GsfElectronToCutBasedSpring15V");
//   TH2F *h2_EGMSFMVA=(TH2F*)f_EGMSF1->Get("GsfElectronToMVAVLooseTightIP2D");
//   TH2F *h2_EGMSFMiniIso=(TH2F*)f_EGMSF1->Get("MVAVLooseElectronToMini");
//   TH2F *h2_EGMSF1=(TH2F*)h2_EGMSFV->Clone("EGMSF1");
//   h2_EGMSF1->Multiply(h2_EGMSFMiniIso);// h2_EGMSF1->Divide(h2_EGMSFMVA);                                                                             
//   TH2F *h2_EGMSF2=(TH2F*)f_EGMSF2->Get("EGamma_SF2D");
//   cout<<"applying EGM SFs to electrons? "<<applyEGMSFs<<endl;


  
  
  double survived_events =0,survived_vetohad=0,not_accepted=0,survived_accept=0,check=0,survived_failed_id=0,survived_failed_iso=0,survived_all=0,events_cr=0,events_id=0,events_failiso=0;
  int event_iso2=0,event_fakephoton=0,event_miniiso=0,event_failminiiso=0,event_failiso2=0,event_failiso=0,event_photonfakes=0,event_iso=0;
  int e_mu_event=0,events_gen_el=0,count_lost_events=0;
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

      /////////////////////////         Weights            /////////////////////////////

      if(Weight<0) continue;  // 2017 and 2018 has negative weights
      Double_t wt =35.9*1000*(Weight);

      
      //in event loop
      double w = puhist->GetBinContent(puhist->GetXaxis()->FindBin(min(TrueNumInteractions,puhist->GetBinLowEdge(puhist->GetNbinsX()+1))));

      wt = wt*w;
      
      // if(strcmp(year,"2016")==0)
      // 	{
      // 	  wt = wt*w;
      // 	}
      /////////////      Overlap removal for Dr(photon,parton)  /////////////////////////////

      // TTGJets and WGJets ---- all events 
            
      if((strcmp(data,"TTJets")==0)&&(madMinPhotonDeltaR > 0.3)) continue;
      if((strcmp(data,"WJets")==0)&&(madMinPhotonDeltaR > 0.5)) continue;

      // TTJets MadHT cut

      // if((strcmp(data,"TTJets")==0)&&(madHT > 600)) continue;

      
	
      

      ////////////////// filters   ///////////////////////////////////////////////
      if(strcmp(year,"2016")==0)
	{
	  if(PrimaryVertexFilter!=1 || globalSuperTightHalo2016Filter!=1 || HBHEIsoNoiseFilter!=1 || HBHENoiseFilter!=1 || EcalDeadCellTriggerPrimitiveFilter!=1)    continue;
	  if(BadPFMuonFilter == false) continue;  
	}
      if(strcmp(year,"2017")==0)
	{
	  if(PrimaryVertexFilter!=1 || globalSuperTightHalo2016Filter!=1 || HBHEIsoNoiseFilter!=1 || HBHENoiseFilter!=1 || EcalDeadCellTriggerPrimitiveFilter!=1)    continue;
	  if(BadPFMuonFilter == false) continue;  
	}

      vector<TLorentzVector> goodphotons;
     
      for(int i=0;i<Photons->size();i++)
	{ if(((*Photons)[i].Pt()>100)&&(TMath::Abs((*Photons)[i].Eta())<2.4)&&((*Photons_fullID)[i]==1)&&((*Photons_hasPixelSeed)[i]<0.000001)/*&&((TMath::Abs((*Photons)[i].Eta())>1.566)&&(TMath::Abs((*Photons)[i].Eta())<1.4446))*/)
	    { goodphotons.push_back((*Photons)[i]);
	      
	    }
	  
	}
      TLorentzVector goodphoton;
      sortTLorVec(&goodphotons);
      if(goodphotons.size()!=0)
	{goodphoton = goodphotons[0];}

      //      if(goodphoton.Pt()<100) continue; // remove events with 0 photons or pt<100
      //if(MET<100) continue;
      
      
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
      //if(dphi1>0.3 && dphi2>0.3) continue;
      //if((st<800 || goodphoton.Pt()<100) && (st<500 || goodphoton.Pt()<190)) continue;
      //if(isoMuonTracks!=0 || isoPionTracks!=0) continue;

      double mt_ele=0,mt_pho=0,mt_elepho=0;
      if(Electrons->size()==1){
	mt_ele=sqrt(2*(*Electrons)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*Electrons)[0].Phi()))));
	if(mt_ele>100) continue;
	if( ((*Electrons)[0].Pt() < 10) || abs((*Electrons)[0].Eta()) > 2.5 ) continue;
      }
      mt_pho=sqrt(2*goodphoton.Pt()*MET*(1-cos(DeltaPhi(METPhi,goodphoton.Phi()))));
      if(Electrons->size()==1) mt_elepho=sqrt(2*(goodphoton+(*Electrons)[0]).Pt()*MET*(1-cos(DeltaPhi(METPhi,(goodphoton+(*Electrons)[0]).Phi()))));   // MHT cut defined in the paper for isolated tracks

    
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

     
      
      ////////         Preselection cut                 //////////////////
      if((goodphoton.Pt()>100)&&(MET>=100)&&(dphi1>0.3)&&(dphi2>0.3)&&(((st>800 && goodphoton.Pt()>100) || (st>500 && goodphoton.Pt()>190)))&&(isoMuonTracks==0)&&(isoPionTracks==0)&&(goodjets.size()>=2)&&(NMuons==0))
	{ 
	  //  cout<<"Not here"<<endl;

	  survived_events+=1;

	  // For prediction using the transfer factors.

	  if(NElectrons == 1)
	    {    fill_hist(cr_pred1,cr_pred2,BTags,MET,goodjets.size(),wt);
	    }  

	   nel[0]->Fill(NElectrons,wt);
	   nmu[0]->Fill(NMuons,wt);
      
	  // Getting the Gen Information

	  // ----------------------  Gen level -------------------------
	  vector<TLorentzVector> gen_el, gen_mu,gen_tau,gen_tau_el,gen_tau_mu,gen_ph,gen_quarks,gen_lep;
     
	  for(int i=0;i<GenParticles->size();i++)
	    { if((*GenParticles)[i].Pt()!=0)
		
		{
		  if((abs((*GenParticles_PdgId)[i])==11) && (abs((*GenParticles_ParentId)[i])==24 || abs((*GenParticles_ParentId)[i])==15) && ((*GenParticles_Status)[i]==1))
		    { 
		      gen_el.push_back((*GenParticles)[i]);
		      if(abs((*GenParticles_ParentId)[i])==15)
			{gen_tau_el.push_back((*GenParticles)[i]);}
		    }
		  if((abs((*GenParticles_PdgId)[i])==13) && (abs((*GenParticles_ParentId)[i])==24 || abs((*GenParticles_ParentId)[i])==15) && ((*GenParticles_Status)[i]==1))
		    { 
		      gen_mu.push_back((*GenParticles)[i]);
		      if(abs((*GenParticles_ParentId)[i])==15)
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
		  if(((abs((*GenParticles_PdgId)[i])>=1)&&(abs((*GenParticles_PdgId)[i])<=6))||(abs((*GenParticles_PdgId)[i])==21))
		    { gen_quarks.push_back((*GenParticles)[i]);}
		  if((abs((*GenParticles_PdgId)[i])==11)||(abs((*GenParticles_PdgId)[i])==13)||(abs((*GenParticles_PdgId)[i])==15))
		    { gen_lep.push_back((*GenParticles)[i]);}
		}
	    }
	  //////////////////////   MIn DR plot to check for the colinear radiation of photon /////////////// 
	  for(int i=0;i<gen_lep.size();i++)
	    { h_mindr_lep_ph->Fill(MinDr(gen_lep[i],*Photons),wt);
	    }
	  for(int i=0;i<Photons->size();i++)
	    { h_mindr_ph_qu->Fill(MinDr((*Photons)[i],gen_quarks),wt);
	    }
	  for(int i=0;i<Photons->size();i++)
	    { h_mindr_ph_lep->Fill(MinDr((*Photons)[i],gen_lep),wt);
	    }
	  for(int i=0;i<gen_quarks.size();i++)
	    { h_mindr_qu_ph->Fill(MinDr(gen_quarks[i],*Photons),wt);
	    }
	  
	  h_madminphotonDR->Fill(madMinPhotonDeltaR,wt);
	  h_mindr_goodph_lep->Fill(MinDr(goodphoton,gen_lep),wt);
	  /////////////////////////  removing the overlap using dr(photon,lepton)   ////////////////

	  if((strcmp(data,"TTJets")==0) && MinDr(goodphoton,gen_lep) > 0.3) continue;
	  if((strcmp(data,"WJets")==0) && MinDr(goodphoton,gen_lep) > 0.5) continue;
	  
	  
	  /////////////////////////////////////////////////////////////////////////////////////////////////
		
	  if(gen_el.size()==0 && gen_mu.size()==0 && gen_tau_el.size()==0 && gen_tau_mu.size()==0) continue; // rejecting all hadronic decays
	  if(gen_el.size()==1 && gen_mu.size()==1)
	    { e_mu_event+=1;}
      
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

      
	  nel[1]->Fill(NElectrons,wt);
	  nmu[1]->Fill(NMuons,wt);
	  //if(gen_el.size()!=1) continue;
	  //if(gen_mu.size()!=0) continue;
	  //if(NElectrons>=2) continue;
	  survived_vetohad+=1;
	  h_lead_ph_pt[1]->Fill(goodphoton.Pt(),wt);
	  h_lead_ph_eta[1]->Fill(goodphoton.Eta(),wt);
  

	  // fill_hist(total1,total2,BTags,MET,goodjets.size(),wt);
      
	  // if(gen_el.size()==1)   // Use to get similar results as of 6 Dec 2019
	  //   {
	  if(gen_mu.size()==0 || (gen_el.size()==1 && gen_mu.size()==1))
	    {
		  
		  
	      // Filling some histograms
	      for(int j=0;j<goodjets.size();j++)
		{ fill_hist_2D(jet_ptbins,BTags,goodjets.size(),goodjets[j].Pt(),wt);
		  fill_hist_2D(mindr2D_el_jet,BTags,goodjets.size(),MinDr(goodjets[j],*Electrons),wt);
		  fill_hist_2D(mindr2D_genel_jet,BTags,goodjets.size(),MinDr(goodjets[j],gen_el),wt);
		      
		}
	      for(int j=0;j<Electrons->size();j++)
		{ fill_hist_2D(e_ptbins,BTags,goodjets.size(),(*Electrons)[j].Pt(),wt);
		}
		  
	      int total_lost_el = 0,cr_el=0;

	      
	      for(int i=0;i<gen_el.size();i++)
		{ 
		  events_gen_el+=1;
		  
		  // Fake photon /////////////////////////////////////////////////
		  double dr=100,match_ph_pt=0,match_el_pt=0;
		  bool isrealphoton = true;
		  bool matched = electron_match_photon(goodphoton);
		  int match_el = 0,match_ph=0;
		  
		  if(gen_el[i].DeltaR(goodphoton)<0.2)
		    { match_el=1;
		      match_el_pt=gen_el[i].Pt();
		    }
		  for(int j=0;j<gen_ph.size();j++)
		    { dr = goodphoton.DeltaR(gen_ph[j]);
		      if(dr<0.2)
			{ match_ph=1;
			  match_ph_pt=gen_ph[j].Pt();
			}
		    }
      
		  if(Electrons->size()==0)
		    { if(match_el==0) isrealphoton=true;
		      else if(match_el==1 && match_ph==0) isrealphoton=false;
		      else if(match_el==1 && match_ph==1)
			{ if(abs(goodphoton.Pt()-match_ph_pt)>abs(goodphoton.Pt()-match_el_pt))
			    { isrealphoton=false;}
			}
		      else isrealphoton=true;
		    }
		  for(int j=0;j<Electrons->size();j++)
		    {  mindr_reco_el_ph[0]->Fill(MinDr((*Electrons)[j],*Photons),wt);}
      
		  mindr_gen_el_reco_ph[0]->Fill(gen_el[i].DeltaR(goodphoton),wt);
		  if(gen_el[i].DeltaR(goodphoton)<0.2)
		    { isrealphoton=false;}
		  if(Photons_electronFakes)
		    {      event_photonfakes+=1;	}
		  if(!isrealphoton || matched)
		    { event_fakephoton+=1;
		      fill_hist(fake_photon1,fake_photon2,BTags,MET,goodjets.size(),wt);
		      //continue;
		    }
		  //if(!isrealphoton) continue;

		  for(int j=0;j<Electrons->size();j++)
		    { mindr_reco_el_ph[1]->Fill(MinDr((*Electrons)[j],*Photons),wt);}
      
		  mindr_gen_el_reco_ph[1]->Fill(gen_el[i].DeltaR(goodphoton),wt);
      
		  
		  

		  fill_hist(total1,total2,BTags,MET,goodjets.size(),wt);
		  bool acceptance = true;
		  bool identified = true;
		  bool isol       = true;
		  bool cr         = true;
		  if(isrealphoton && !matched)
		    {
      
		      // Out of acceptance ////////////////////////////////////////////////
		      
		      acceptance = true;
		      if((abs(gen_el[i].Eta())>2.5)||(gen_el[i].Pt()<10)&&(gen_mu.size()!=1)||(Electrons->size()==0))
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
			  identified = true;
			  mindr_gen_rec_el ->Fill(MinDr(gen_el[i],*Electrons),wt);
			  
			  if((MinDr(gen_el[i],*Electrons)>0.2)&&(gen_mu.size()!=1))/* || (gen_el.size()>0 && Electrons->size()==0)))*/ 
			    { events_id+=1;
			      fill_hist(fail_id1,fail_id2,BTags,MET,goodjets.size(),wt);
			      identified = false;
			    }
			}
		      if(identified && acceptance) 
			{
			  ////////////////////////////////////////////////////////////
			  survived_failed_id+=1;
			  nel[3]->Fill(NElectrons,wt);
			  nmu[3]->Fill(NMuons,wt);
			  gen_el_size[2]->Fill(gen_el.size(),wt);
			  gen_mu_size[2]->Fill(gen_mu.size(),wt);
			  gen_ph_size[2]->Fill(gen_ph.size(),wt);
			  
			      
			  
			  //   Failed Iso ///////////////////////////////////////////
			  
			  isol= true;
			  if((NElectrons == 0) &&(gen_mu.size()!=1))
			    { event_failiso+=1; 
			      fill_hist(fail_iso1,fail_iso2,BTags,MET,goodjets.size(),wt);
			      isol = false;
			    }
			  if(Electrons->size()!=0 && !Electrons_passIso)
			    { event_failiso2+=1;	}
			  if(NElectrons>0)
			    { event_iso+=1;} 
			  if(Electrons->size()!=0 && Electrons_passIso)
			    { event_iso2+=1;	}
			  
			  
			  bool miniIso=false;
			  for(int j=0;j<Electrons_MiniIso->size();j++)
			    { if((*Electrons_MiniIso)[j]<0.1)
				{ miniIso = true;}
			    }
			  if(miniIso)
			    { event_miniiso+=1;}
			  else if(!miniIso)
			    { event_failminiiso+=1;}
			}
		      if(!acceptance || !isol || !identified)
			{
			  total_lost_el += 1;  
			}
		      
		      if(isol && acceptance && identified) 
		      	{
			   /////////////////////////////////////////////////////////
			  survived_failed_iso+=1;
			  nel[4]->Fill(NElectrons,wt);
			  nmu[4]->Fill(NMuons,wt);
			  gen_el_size[3]->Fill(gen_el.size(),wt);
			  gen_mu_size[3]->Fill(gen_mu.size(),wt);
			  gen_ph_size[3]->Fill(gen_ph.size(),wt);
			}
		      /////////   } // is real photon
		    	  
			  // 1 Lepton CR ///////////////////////////////////////////
			  cr = true;
			  if((NElectrons == 1)&&(cr_el<1))
			    { events_cr+=1;
			      fill_hist(one_lep_cr1,one_lep_cr2,BTags,MET,goodjets.size(),wt);
			      //	      SR_hist(srbins_cr,BTags,MET,goodjets.size(),wt,year,1);
			      if(BTags==0)
				{
				  METNJ_B0_E1->Fill(MET,goodjets.size(),wt);
				}
			      if(BTags>=1)
				{
				  METNJ_B1_E1->Fill(MET,goodjets.size(),wt);
				}
			      cr = false;
			      cr_el+=1;
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
			  int njets = goodjets.size();
			  int btags = BTags;
			  if(btags==0)
			    {
			      if(njets>=2 && njets<=4)
				{ for(int k=0;k<goodjets.size();k++)
				    {
				      jetpt_inbins[0]->Fill(goodjets[k].Pt(),wt);
				      mindr1D_el_jet[0]->Fill(MinDr(goodjets[k],*Electrons),wt);
				      mindr1D_genel_jet[0]->Fill(MinDr(goodjets[k],gen_el),wt);
				    }
				  for(int k=0;k<Electrons->size();k++)
				    {
				      ept_inbins[0]->Fill((*Electrons)[k].Pt(),wt);
				    }
				}
			      if(njets>=5 && njets<=6)
				{
				  for(int k=0;k<goodjets.size();k++)
				    {
				      jetpt_inbins[1]->Fill(goodjets[k].Pt(),wt);
				      mindr1D_el_jet[1]->Fill(MinDr(goodjets[k],*Electrons),wt);
				      mindr1D_genel_jet[1]->Fill(MinDr(goodjets[k],gen_el),wt);
				    }
				  for(int k=0;k<Electrons->size();k++)
				    {
				      ept_inbins[1]->Fill((*Electrons)[k].Pt(),wt);
				    }
				}
			      if(njets>=7)
				{
				  for(int k=0;k<goodjets.size();k++)
				    {
				      jetpt_inbins[2]->Fill(goodjets[k].Pt(),wt);
				      mindr1D_el_jet[2]->Fill(MinDr(goodjets[k],*Electrons),wt);
				      mindr1D_genel_jet[2]->Fill(MinDr(goodjets[k],gen_el),wt);
				    }
				  for(int k=0;k<Electrons->size();k++)
				    {
				      ept_inbins[2]->Fill((*Electrons)[k].Pt(),wt);
				    }
				}

			    }
			  if(btags==1)
			    { if(njets>=2 && njets<=4)
				{
				  for(int k=0;k<goodjets.size();k++)
				    {
				      jetpt_inbins[3]->Fill(goodjets[k].Pt(),wt);
				      mindr1D_el_jet[3]->Fill(MinDr(goodjets[k],*Electrons),wt);
				      mindr1D_genel_jet[3]->Fill(MinDr(goodjets[k],gen_el),wt);
				    }
				  for(int k=0;k<Electrons->size();k++)
				    {
				      ept_inbins[3]->Fill((*Electrons)[k].Pt(),wt);
				    }
				}
			      if(njets>=5 && njets<=6)
				{for(int k=0;k<goodjets.size();k++)
				    {
				      jetpt_inbins[4]->Fill(goodjets[k].Pt(),wt);
				      mindr1D_el_jet[4]->Fill(MinDr(goodjets[k],*Electrons),wt);
				      mindr1D_genel_jet[4]->Fill(MinDr(goodjets[k],gen_el),wt);
				    }
				  for(int k=0;k<Electrons->size();k++)
				    {
				      ept_inbins[4]->Fill((*Electrons)[k].Pt(),wt);
				    }
					     
				}
			      if(njets>=7)
				{
				  for(int k=0;k<goodjets.size();k++)
				    {
				      jetpt_inbins[5]->Fill(goodjets[k].Pt(),wt);
				      mindr1D_el_jet[5]->Fill(MinDr(goodjets[k],*Electrons),wt);
				      mindr1D_genel_jet[5]->Fill(MinDr(goodjets[k],gen_el),wt);
				    }
				  for(int k=0;k<Electrons->size();k++)
				    {
				      ept_inbins[5]->Fill((*Electrons)[k].Pt(),wt);
				    }
				}
	      
			    }
			  if(btags>=2)
			    { if(njets>=2 && njets<=4)
				{
				  for(int k=0;k<goodjets.size();k++)
				    {
				      jetpt_inbins[6]->Fill(goodjets[k].Pt(),wt);
				      mindr1D_el_jet[6]->Fill(MinDr(goodjets[k],*Electrons),wt);
				      mindr1D_genel_jet[6]->Fill(MinDr(goodjets[k],gen_el),wt);
				    }
				  for(int k=0;k<Electrons->size();k++)
				    {
				      ept_inbins[6]->Fill((*Electrons)[k].Pt(),wt);
				    }
				}
			      if(njets>=5 && njets<=6)
				{
				  for(int k=0;k<goodjets.size();k++)
				    {
				      jetpt_inbins[7]->Fill(goodjets[k].Pt(),wt);
				      mindr1D_el_jet[7]->Fill(MinDr(goodjets[k],*Electrons),wt);
				      mindr1D_genel_jet[7]->Fill(MinDr(goodjets[k],gen_el),wt);
				    }
				  for(int k=0;k<Electrons->size();k++)
				    {
				      ept_inbins[7]->Fill((*Electrons)[k].Pt(),wt);
				    }
				}
			      if(njets>=7)
				{
				  for(int k=0;k<goodjets.size();k++)
				    {
				      jetpt_inbins[8]->Fill(goodjets[k].Pt(),wt);
				      mindr1D_el_jet[8]->Fill(MinDr(goodjets[k],*Electrons),wt);
				      mindr1D_genel_jet[8]->Fill(MinDr(goodjets[k],gen_el),wt);
				    }
				  for(int k=0;k<Electrons->size();k++)
				    {
				      ept_inbins[8]->Fill((*Electrons)[k].Pt(),wt);
				    }
				}
	      
			    }
					  
			} // cr
		      /*	} //fail iso
				} // fail id
				} // fail accept
		      */	} //real photon
		
	    
		 
		}// loop over gen_electrons
	      if(cr_el>=2)
	      	{ cout<<"problem"<<endl;}
	      if(total_lost_el==2 || (total_lost_el==1)/* && cr_el!=1) */)
	      	{ 
	      	  fill_hist(lost_event1,lost_event2,BTags,MET,goodjets.size(),wt);
	      	  count_lost_events +=1;
		  //		  SR_hist(srbins_lostel,BTags,MET,goodjets.size(),wt,year,0);
		  if(BTags==0)
		    {  METNJ_B0_E0->Fill(MET,goodjets.size(),wt);}
		  if(BTags>=1)
		    {  METNJ_B1_E0->Fill(MET,goodjets.size(),wt);}
		}
	    
	      //}// end of gen_mu condition
	    }// end of gen_el.size()!=0 condition
	 
	}// end of preselection cuts
    }// loop over entries
  
  cout<<"Events survived preselection             = "<<survived_events<<endl;
  cout<<"Events survived veto had                 = "<<survived_vetohad<<endl;
  cout<<"Events with gen e and gen mu after presel= "<<e_mu_event<<endl;
  cout<<"Events with atleast 1 gen electron       = "<<events_gen_el<<endl;
  cout<<"Events where gen electron fakes as photon= "<<event_fakephoton<<endl;
  cout<<"Events passing Photon_electronfakes      = "<<event_photonfakes<<endl;
  cout<<"Events not accepted by the detector      = "<<not_accepted<<endl;
  cout<<"Events survived acceptance cut           = "<<survived_accept<<endl;
  cout<<"Events failed ID                         = "<<events_id<<endl;
  cout<<"Events survived id cut                   = "<<survived_failed_id<<endl;
  cout<<"Events failed Iso (NElecltrons==0)       = "<<event_failiso<<endl;
  cout<<"Events failed pass iso (!passiso)        = "<<event_failiso2<<endl;
  cout<<"Events failed mini iso (MiniIso>0.2)     = "<<event_failminiiso<<endl;
  cout<<"Events pass NElectrons==0                = "<<event_iso<<endl;
  cout<<"Events pass (passIso)                    = "<<event_iso2<<endl;
  cout<<"Events pass mini_iso <0.1                = "<<event_miniiso<<endl;
  cout<<"Events survived iso cut                  = "<<survived_failed_iso<<endl;
  cout<<"Events in 1 lepton CR                    = "<<events_cr<<endl;
  cout<<"Events survived all cut                  = "<<survived_all<<endl;
  cout<<endl;
  cout<<"Events which are lost                    = "<<count_lost_events<<endl;

  // TF_B0->Copy(
  // TF_B0 = METNJ_B0_E0->Divide(METNJ_B0_E1);
  // TF_B1 = METNJ_B1_E1->Divide(METNJ_B1_E1);
}
	
  
bool lost_el::electron_match_photon(TLorentzVector photon)
{ for(int i=0;i<Electrons->size();i++)
    { if(photon.DeltaR((*Electrons)[i])<0.2)
	{ return true;	}
      else
	{ return false;}
    }
}

void lost_el::fill_hist(TH1D *h1, TH1D *h2,int btags,double met,int njets,double wt)
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


void lost_el::fill_hist_2D(TH2D *h1,int btags,int njets,double v1,double wt)
{
  //h1 = TH1D("h1","lost electron in b-jets and njets bins",6,1,7);
  //h2 = TH1D("h2","lost electron in all the bins",16,1,17);
  
  // h1->Fill("NJets_{=0}^{2-4}",0,0);
  // h1->Fill("NJets_{= 1}^{2-4}",0);
  // h1->Fill("NJets_{#geq 2}^{2-4}",0);
  // h1->Fill("NJets_{=0}^{5-6}",0);
  // h1->Fill("NJets_{= 1}^{5-6}",0);
  // h1->Fill("NJets_{#geq 2}^{5-6}",0);
  // h1->Fill("NJets_{=0}^{#geq 7}",0);
  // h1->Fill("NJets_{= 1}^{#geq 7}",0);
  // h1->Fill("NJets_{#geq 2}^{#geq 7}",0);
  
  //int njets = goodjets.size();
  
  if(btags==0)
    {
      if(njets>=2 && njets<=4)
	{		  
	  h1->Fill("NJets_{=0}^{2-4}",v1,wt);
	}
      if(njets>=5 && njets<=6)
	{
	  h1->Fill("NJets_{=0}^{5-6}",v1,wt);
	}
      if(njets>=7)
	{
	  h1->Fill("NJets_{=0}^{#geq 7}",v1,wt);
	}

    }
  if(btags==1)
    { if(njets>=2 && njets<=4)
	{
	  h1->Fill("NJets_{= 1}^{2-4}",v1,wt);
	}
      if(njets>=5 && njets<=6)
	{
	  h1->Fill("NJets_{= 1}^{5-6}",v1,wt);
	}
      if(njets>=7)
	{
	  h1->Fill("NJets_{= 1}^{#geq 7}",v1,wt);
	}
	      
    }
  if(btags>=2)
    { if(njets>=2 && njets<=4)
	{
	  h1->Fill("NJets_{#geq 2}^{2-4}",v1,wt);
	}
      if(njets>=5 && njets<=6)
	{
	  h1->Fill("NJets_{#geq 2}^{5-6}",v1,wt);
	}
      if(njets>=7)
	{
	  h1->Fill("NJets_{#geq 2}^{#geq 7}",v1,wt);
	}
	      
    }
}

// void lost_el::SR_hist(TH1D *h1,int btags,double met,int njets,double wt, const char *year, int x)
// {
//   //h1 = TH1D("h1","lost electron in b-jets and njets bins",6,1,7);
//   //h2 = TH1D("h2","lost electron in all the bins",16,1,17);
//   TFile *f1;
//   if(strcmp(year,"2016")==0)
//     {
//       f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2016_Transfer_Factors.root","READ");
//     }
//   if(strcmp(year,"2017")==0)
//     {
//       f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2017_Transfer_Factors.root","READ");
//     }
//   if(strcmp(year,"2018")==0)
//     {
//       f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2018_Transfer_Factors.root","READ");
//     }
//   TH1D *tf;
//   tf = (TH1D*)f1->Get("tf");
//   if(x=1)
//     { if(btags==0)
// 	{ if(njets == 2){ wt=wt*tf->GetBinContent(2);}
// 	  if(njets == 3){ wt=wt*tf->GetBinContent(4);}
// 	  if(njets == 4){ wt=wt*tf->GetBinContent(6);}
// 	  if(njets>=5 && njets<=6){ wt=wt*tf->GetBinContent(8);}
// 	  if(njets >= 7){ wt=wt*tf->GetBinContent(10);}
// 	}
//       if(btags>=1)
// 	{ if(njets>=2 && njets<=4){ wt=wt*tf->GetBinContent(12);}
// 	  if(njets>=5 && njets<=6){ wt=wt*tf->GetBinContent(14);}
// 	  if(njets>=7){ wt=wt*tf->GetBinContent(16);}
// 	}
//     }
	  
	  
  
//   if(btags==0)
//     { if(njets>= 2 && njets<=4)
// 	{ if(met>=200 && met <270)
// 	    { h1->Fill(1,wt);}
// 	  if(met>=270 && met <350)
// 	    { h1->Fill(2,wt);}
// 	  if(met>=350 && met <450)
// 	    { h1->Fill(3,wt);}
// 	  if(met>=450 && met <750)
// 	    { h1->Fill(4,wt);}
// 	  if(met>=750)
// 	    { h1->Fill(5,wt);}
// 	}
//       if(njets>= 5 && njets<=6)
// 	{ if(met>=200 && met <270)
// 	    { h1->Fill(6,wt);}
// 	  if(met>=270 && met <350)
// 	    { h1->Fill(7,wt);}
// 	  if(met>=350 && met <450)
// 	    { h1->Fill(8,wt);}
// 	  if(met>=450)
// 	    { h1->Fill(9,wt);}
// 	}
//       if(njets>= 7)
// 	{ if(met>=200 && met <270)
// 	    { h1->Fill(10,wt);}
// 	  if(met>=270 && met <350)
// 	    { h1->Fill(11,wt);}
// 	  if(met>=350 && met <450)
// 	    { h1->Fill(12,wt);}
// 	  if(met>=450)
// 	    { h1->Fill(13,wt);}
// 	}
//     }
//   if(btags>=1)
//     { if(njets>= 2 && njets<=4)
// 	{ if(met>=200 && met <270)
// 	    { h1->Fill(14,wt);}
// 	  if(met>=270 && met <350)
// 	    { h1->Fill(15,wt);}
// 	  if(met>=350 && met <450)
// 	    { h1->Fill(16,wt);}
// 	  if(met>=450)
// 	    { h1->Fill(17,wt);}
// 	}
//       if(njets>= 5 && njets<=6)
// 	{ if(met>=200 && met <270)
// 	    { h1->Fill(18,wt);}
// 	  if(met>=270 && met <350)
// 	    { h1->Fill(19,wt);}
// 	  if(met>=350 && met <450)
// 	    { h1->Fill(20,wt);}
// 	  if(met>=450)
// 	    { h1->Fill(21,wt);}
// 	}
//       if(njets>= 7)
// 	{ if(met>=200 && met <270)
// 	    { h1->Fill(22,wt);}
// 	  if(met>=270 && met <350)
// 	    { h1->Fill(23,wt);}
// 	  if(met>=350 && met <450)
// 	    { h1->Fill(24,wt);}
// 	  if(met>=450)
// 	    { h1->Fill(25,wt);}
// 	}
//     }
// }
