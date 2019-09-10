void TF_lost_el_all_region(char* input)
{
  TFile *f1,*f2,*f3,*f4,*f5,*f6;

  
  if(strcmp(input,"wgjets_lnu")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_ll_mu.root");
    }
  else if(strcmp(input,"ttjets")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/TTJets/ttjets_ll_mu.root");
    }

  else if(strcmp(input,"ttgjets")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/TTGJets/ttgjets_ll_mu.root");
    }

  else if(strcmp(input,"gjets")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/GJets/gjets_ll_mu.root");
    }

  else if(strcmp(input,"wgjets_monophoton")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_ll_mu.root");
    }
  else if(strcmp(input,"wgjets")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/WGJets/wgjets_ll_mu.root");
    }
  else
    { cout<<"please give a proper input"; return 0;
    }


  TH1D *fail_accept,*fail_id,*fail_iso,*cr,*fake_photon,*total;
  
  fake_photon = (TH1D*)f1->Get("fake_photon_2");
  fail_accept = (TH1D*)f1->Get("fail_accept_2");
  fail_id     = (TH1D*)f1->Get("fail_id_2");
  fail_iso    = (TH1D*)f1->Get("fail_iso_2");
  cr          = (TH1D*)f1->Get("one_lep_cr_2");
  total       = (TH1D*)f1->Get("total2");
  
  
  
  
  // double total=0;
  const char* str[] = {"NJets_{0}^{=2} & 100<MET<150","NJets_{0}^{=2} & MET#geq 150","NJets_{0}^{=3} & 100<MET<150","NJets_{0}^{=3} & MET#geq 150","NJets_{0}^{=4} & 100<MET<150","NJets_{0}^{=4} & MET#geq 150","NJets_{0}^{5-6} & 100<MET<150","NJets_{0}^{5-6} & MET#geq 150","NJets_{0}^{#geq7} & 100<MET<150","NJets_{0}^{#geq7} & MET#geq 150","NJets_{#geq 1}^{2-4} & 100<MET<150","NJets_{#geq 1}^{2-4} & MET#geq 150","NJets_{#geq 1}^{5-6} & 100<MET<150","NJets_{#geq 1}^{5-6} & MET#geq 150","NJets_{#geq 1}^{#geq 7} & 100<MET<150","NJets_{#geq 1}^{#geq 7} & MET#geq 150"};
  
  // for(int i=1;i<=6;i++)
  //   {
  //     total=fake_photon->GetBinContent(i)+had_tau->GetBinContent(i)+fail_accept->GetBinContent(i)+fail_id->GetBinContent(i)+fail_iso->GetBinContent(i)+cr->GetBinContent(i);
  //     fake_photon->SetBinContent(i,fake_photon->GetBinContent(i)/total);
  //     had_tau->SetBinContent(i,had_tau->GetBinContent(i)/total);
  //     fail_accept->SetBinContent(i,fail_accept->GetBinContent(i)/total);
  //     fail_id->SetBinContent(i,fail_id->GetBinContent(i)/total);
  //     fail_iso->SetBinContent(i,fail_iso->GetBinContent(i)/total);
  //     cr->SetBinContent(i,cr->GetBinContent(i)/total);
  //   }


  // Making new copy of hist with total divided

  fake_photon->Divide(total);
  fail_accept->Divide(total);
  fail_id->Divide(total);
  fail_iso->Divide(total);
  cr->Divide(total);
  
  THStack *stack = new THStack("Stack","stack hist");
  TCanvas *c1 = new TCanvas("stackhist","stackhist",1600,900);

  
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);

  TPad *pad2 = new TPad("pad1","pad1",0,0.0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad1->Draw();pad1->SetGridx();
  pad2->Draw();pad2->SetGridx();pad2->SetGridy();
  pad1->cd();
  
  gStyle->SetPalette(kOcean);
  fail_accept->SetFillStyle(3008);
  fail_accept->SetFillColor(kGreen);
  fail_id->SetFillStyle(3019);
  fail_id->SetFillColor(kOrange);
  fail_iso->SetFillStyle(3144);
  fail_iso->SetFillColor(kBlue+3);
  fake_photon->SetFillColor(kBlue-9);
  cr->SetFillColor(kGray+3);
  
  TLegend *legend = new TLegend(0.7,0.2,0.9,0.4);
  stack->Add(cr);
  stack->Add(fail_accept);
  stack->Add(fail_id);
  stack->Add(fail_iso);
  stack->Add(fake_photon);
  stack->Draw("hist");
  legend->SetNColumns(3);
  legend->SetBorderSize(1);

  stack->GetYaxis()->SetTitleSize(0.05);
  stack->GetYaxis()->SetTitleOffset(0.85);
  stack->GetYaxis()->SetLabelSize(0.07);
  stack->GetYaxis()->SetTitle(0);
  stack->GetYaxis()->SetRangeUser(0,2);
  
  stack->SetTitle(input);


  legend->AddEntry(cr,"1L CR","f");
  legend->AddEntry(fail_accept,"Fail Accept","f");
  legend->AddEntry(fail_id,"Fail Id","f");
  legend->AddEntry(fail_iso,"Fail Iso","f");
  legend->AddEntry(fake_photon,"Fake photon","f");
  legend->SetTextSize(0.03);
  legend->Draw();


  pad2->cd();

  TH1D *TF = new TH1D("tf","Transfer factor",16,1,17);
  for(int i=1;i<17;i++)
    { TF->GetXaxis()->SetBinLabel(i,str[i-1]); }
  TF->Add(fail_accept);
  TF->Add(fail_id);
  TF->Add(fail_iso);
  //TF->Add(fake_photon_1);
  TF->GetYaxis()->SetRangeUser(0,2);
  TF->Sumw2();
  TF->SetStats(0);
  TF->Divide(cr);
  TF->Draw("ep");
  TF->SetTitle(0);

  TF->GetXaxis()->SetTitle(0);
  TF->GetXaxis()->SetLabelSize(0.06);

  TF->GetYaxis()->SetTitle("Transfer factor muons");
  TF->GetYaxis()->SetTitleOffset(0.35);
  TF->GetYaxis()->SetTitleSize(0.13);
  TF->GetYaxis()->SetLabelSize(0.09);
  TF->SetLineWidth(3);
  //TF->SetMaximum(1.99);
  //TF->SetMinimum(0.01);

  for(int i=1;i<17;i++)
    { cout<<"The Transfer Factor in bin"<<i<<" = "<<TF->GetBinContent(i)<<" +- "<<TF->GetBinError(i)<<endl; }
}

