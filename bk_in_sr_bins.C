void bk_in_sr_bins(char* input)
{
  TFile *f1,*f2,*f3,*f4,*f5,*f6;
  if(strcmp(input,"wgjets_lnu")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_fake_photon.root");
      f2 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_had_tau.root");
      f3 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_ll_mu.root");
      f4 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_ll.root");
      f5 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_one_el.root");
      f6 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_one_mu.root");
    }
  else if(strcmp(input,"ttjets")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/TTJets/ttjets_fake_photon.root");
      f2 = new TFile("../../EHEP_Projects/Summer_2019/TTJets/ttjets_had_tau.root");
      f3 = new TFile("../../EHEP_Projects/Summer_2019/TTJets/ttjets_ll_mu.root");
      f4 = new TFile("../../EHEP_Projects/Summer_2019/TTJets/ttjets_ll.root");
      f5 = new TFile("../../EHEP_Projects/Summer_2019/TTJets/ttjets_one_el.root");
      f6 = new TFile("../../EHEP_Projects/Summer_2019/TTJets/ttjets_one_mu.root");
    }

  else if(strcmp(input,"ttgjets")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/TTGJets/ttgjets_fake_photon.root");
      f2 = new TFile("../../EHEP_Projects/Summer_2019/TTGJets/ttgjets_had_tau.root");
      f3 = new TFile("../../EHEP_Projects/Summer_2019/TTGJets/ttgjets_ll_mu.root");
      f4 = new TFile("../../EHEP_Projects/Summer_2019/TTGJets/ttgjets_ll.root");
      f5 = new TFile("../../EHEP_Projects/Summer_2019/TTGJets/ttgjets_one_el.root");
      f6 = new TFile("../../EHEP_Projects/Summer_2019/TTGJets/ttgjets_one_mu.root");
    }

  else if(strcmp(input,"gjets")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/GJets/gjets_fake_photon.root");
      f2 = new TFile("../../EHEP_Projects/Summer_2019/GJets/gjets_had_tau.root");
      f3 = new TFile("../../EHEP_Projects/Summer_2019/GJets/gjets_ll_mu.root");
      f4 = new TFile("../../EHEP_Projects/Summer_2019/GJets/gjets_ll.root");
      f5 = new TFile("../../EHEP_Projects/Summer_2019/GJets/gjets_one_el.root");
      f6 = new TFile("../../EHEP_Projects/Summer_2019/GJets/gjets_one_mu.root");
    }

  else if(strcmp(input,"wgjets_monophoton")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_fake_photon.root");
      f2 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_had_tau.root");
      f3 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_ll_mu.root");
      f4 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_ll.root");
      f5 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_one_el.root");
      f6 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_one_mu.root");

    }

  else
    { cout<<"please give a proper input";
    }
  
  TH1D *ll_e,*ll_mu,*fake_photon,*had_tau,*one_el,*one_mu;

  TH1D *ll,*fake_photon_1,*had_tau_1,*one_lep;
  ll = new TH1D("ll","lost lepton",6,1,7);
  fake_photon_1 = new TH1D("fake_photon_1","fake photon",6,1,7);
  had_tau_1 = new TH1D("had_tau_1","hadronic tau",6,1,7);
  one_lep = new TH1D("one_lep","1 lepton",6,1,7);
  
  fake_photon = (TH1D*)f1->Get("faking_photon");
  had_tau = (TH1D*)f2->Get("had_tau");
  ll_mu = (TH1D*)f3->Get("ll_mu");
  ll_e = (TH1D*)f4->Get("ll_electron");
  one_el = (TH1D*)f5->Get("one_el");
  one_mu = (TH1D*)f6->Get("one_mu");

  // TH1D *total;
  double temp=0,temp_one_lep=0,temp_lost_lep=0;
  const char* str[6] = {"NJets_{=0}^{2-4}","NJets_{#geq 1}^{2-4}","NJets_{=0}^{5-6}","NJets_{#geq 1}^{5-6}","NJets_{=0}^{#geq 7}","NJets_{#geq 1}^{#geq 7}"};
  
  for(int i=1;i<=6;i++)
    { temp_one_lep = one_el->GetBinContent(i)+one_mu->GetBinContent(i);
      temp_lost_lep = ll_mu->GetBinContent(i) + ll_e->GetBinContent(i);
      temp = fake_photon->GetBinContent(i)+had_tau->GetBinContent(i)+temp_one_lep+temp_lost_lep;
      cout<<"temp = "<<temp<<" temp_one_lep = "<<temp_one_lep<<" temp_lost_lep = "<<temp_lost_lep<<endl;
      ll->Fill(str[i-1],temp_lost_lep/temp);
      fake_photon_1->Fill(str[i-1],fake_photon->GetBinContent(i)/temp);
      had_tau_1->Fill(str[i-1],had_tau->GetBinContent(i)/temp);
      one_lep->Fill(str[i-1],temp_one_lep/temp);
    }

  THStack *stack = new THStack("Stack","stack hist");

  TCanvas *c1 = new TCanvas("stackhist","stackhist",1600,900);
  
  gStyle->SetPalette(kOcean);

  ll->SetFillColor(kOrange+3);
  fake_photon_1->SetFillColor(kBlue-9);
  had_tau_1->SetFillColor(kGray);
  one_lep->SetFillColor(kYellow);

  TLegend *legend = new TLegend(0.55,0.7,0.9,0.9);

  stack->Add(one_lep);
  stack->Add(ll);
  stack->Add(had_tau_1);
  stack->Add(fake_photon_1);
  

  stack->Draw("hist");

  legend->SetNColumns(2);
  legend->SetBorderSize(0);

  stack->GetXaxis()->SetTitleSize(0.045);
  stack->GetXaxis()->SetLabelSize(0.045);
  stack->GetXaxis()->SetTitleOffset(0.88);
  stack->GetXaxis()->SetTitle("Bin Number");
  //stack->GetXaxis()->SetGridStyle(1);
  
  stack->GetYaxis()->SetTitleSize(0.05);
  stack->GetYaxis()->SetTitleOffset(0.85);
  stack->GetYaxis()->SetLabelSize(0.05);
  stack->GetYaxis()->SetTitle(0);
  stack->GetYaxis()->SetRange(0,2);
  stack->SetTitle(input);


  legend->AddEntry(one_lep,"1L CR","f");
  legend->AddEntry(ll,"lost lepton","f");
  legend->AddEntry(fake_photon_1,"Fake photon","f");
  legend->AddEntry(had_tau_1,"hadronic #tau","f");
  legend->SetTextSize(0.02);
  legend->Draw();

  TArrow *arrow1 = new TArrow(0.5,1300,6.5,1300,0.01,"<|>");
  arrow1->Draw("");

  TLatex T1;
  T1.SetTextSize(0.08);
  T1.DrawLatex(3.0,2000,"N^{0}_{2-4}");
  
}
