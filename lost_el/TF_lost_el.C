void TF_lost_el(char* input )
{
  TFile *f1,*f2,*f3,*f4,*f5,*f6;

  
  if(strcmp(input,"wgjets_lnu")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_ll.root");
    }
  else if(strcmp(input,"ttjets")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/TTJets/ttjets_ll.root");
    }

  else if(strcmp(input,"ttgjets")==0)
    {
      f1 = new TFile("output_files/TTGJets_lost_el.root");
    }

  else if(strcmp(input,"gjets")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/GJets/gjets_ll.root");
    }

  else if(strcmp(input,"wgjets_monophoton")==0)
    {
      f1 = new TFile("../../EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_ll.root");
    }
  else if(strcmp(input,"wgjets")==0)
    {
      f1 = new TFile("output_files/WGJets_lost_el.root");
    }
  else
    { cout<<"please give a proper input"; return 0;
    }
  
  const char* str[9] = {"NJets_{=0}^{2-4}","NJets_{= 1}^{2-4}","NJets_{#geq 2}^{2-4}","NJets_{=0}^{5-6}","NJets_{= 1}^{5-6}","NJets_{#geq 2}^{5-6}","NJets_{=0}^{#geq 7}","NJets_{= 1}^{#geq 7}","NJets_{#geq 2}^{#geq 7}"};
  

  TH1D *fail_accept,*fail_id,*fail_iso,*cr,*fake_photon;
  
  fake_photon = (TH1D*)f1->Get("fake_photon_1");
  fail_accept = (TH1D*)f1->Get("fail_accept_1");
  fail_id     = (TH1D*)f1->Get("fail_id_1");
  fail_iso    = (TH1D*)f1->Get("fail_iso_1");
  cr          = (TH1D*)f1->Get("one_lep_cr_1");
  TH1D *total = new TH1D("total","Total = fail_id+fail_iso+fail_accept+1e_cr",9,1,10);
  for(int i=1;i<=9;i++)
    { total->GetXaxis()->SetBinLabel(i,str[i-1]);}
  total->Add(fail_accept);
  total->Add(fail_id);
  total->Add(fail_iso);
  total->Add(cr);


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
  fail_id->SetFillStyle(3021);
  fail_id->SetFillColor(kRed);
  fail_iso->SetFillStyle(3144);
  fail_iso->SetFillColor(kBlue+3);
  fake_photon->SetFillColor(kBlue-9);
  //  cr->SetFillColor(kGray);
  
  TLegend *legend = new TLegend(0.7,0.2,0.9,0.4);
  stack->Add(cr);
  stack->Add(fail_accept);
  stack->Add(fail_id);
  stack->Add(fail_iso);
  // stack->Add(fake_photon);
  stack->Draw("hist");
  legend->SetNColumns(2);
  legend->SetBorderSize(1);

  stack->GetYaxis()->SetTitleSize(0.05);
  stack->GetYaxis()->SetTitleOffset(0.85);
  stack->GetYaxis()->SetLabelSize(0.07);
  stack->GetYaxis()->SetTitle(0);
  stack->GetYaxis()->SetRangeUser(0,2);
  char *title = new char[200];
  sprintf(title,"%s_MET_200",input);
  stack->SetTitle(title);


  legend->AddEntry(cr,"1e CR","f");
  legend->AddEntry(fail_accept,"Fail Accept","f");
  legend->AddEntry(fail_id,"Fail Id","f");
  legend->AddEntry(fail_iso,"Fail Iso","f");
  //legend->AddEntry(fake_photon,"Fake photon","f");
  legend->SetTextSize(0.03);
  legend->Draw();


  pad2->cd();

  TH1D *TF = new TH1D("tf","Transfer factor",9,1,10);
  for(int i=1;i<=9;i++)
    { TF->GetXaxis()->SetBinLabel(i,str[i-1]);}
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
  TF->GetXaxis()->SetLabelSize(0.20);

  TF->GetYaxis()->SetTitle("Transfer factor");
  TF->GetYaxis()->SetTitleOffset(0.35);
  TF->GetYaxis()->SetTitleSize(0.13);
  TF->GetYaxis()->SetLabelSize(0.09);
  TF->SetLineWidth(3);
  //TF->SetMaximum(1.99);
  //TF->SetMinimum(0.01);

  for(int i=1;i<=9;i++)
    { cout<<"The Transfer Factor in bin"<<i<<" = "<<TF->GetBinContent(i)<<" +- "<<TF->GetBinError(i)<<endl;}
  
}

