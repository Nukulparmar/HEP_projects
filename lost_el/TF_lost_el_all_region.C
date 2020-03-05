void TF_lost_el_all_region(char* input, int year,char* save,char* save2)
{
  TFile *f1/*,*f2*/;
  // char *temp = new char[200] ;
  if(year == 2016)
    { 
      if(strcmp(input,"ttgjets")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el.root");
	  //f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2016_Transfer_Factors.root","recreate");
	}
      else if(strcmp(input,"wgjets")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2016_Transfer_Factors.root","recreate"); 
	}
      else if(strcmp(input,"both")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTWGamma_lost_el_2016.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTWGamma_lost_el_2016_Transfer_Factors.root","recreate");
	}
      else if(strcmp(input,"all")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2016.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2016_Transfer_Factors.root","recreate");
	}
      else
	{ cout<<"please give a proper input"; return 0;
	}
    }
  if(year == 2017)
    { 
      if(strcmp(input,"ttgjets")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2017.root");
	  //  f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2017_Transfer_Factors.root","recreate");
	}
      else if(strcmp(input,"wgjets")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2017.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2017_Transfer_Factors.root","recreate");
	}
      else if(strcmp(input,"both")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTWGamma_lost_el_2017.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTWGamma_lost_el_2017_Transfer_Factors.root","recreate");
	}
      else if(strcmp(input,"all")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2017.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2017_Transfer_Factors.root","recreate");
	}
      else
	{ cout<<"please give a proper input"; return 0;
	}
    }
  if(year == 2018)
    { 
      if(strcmp(input,"ttgjets")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2018.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2018_Transfer_Factors.root","recreate");
	}
      else if(strcmp(input,"wgjets")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2018.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2018_Transfer_Factors.root","recreate");
	}
      else if(strcmp(input,"both")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTWGamma_lost_el_2018.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/TTWGamma_lost_el_2018_Transfer_Factors.root","recreate");
	}
      else if(strcmp(input,"all")==0)
	{
	  f1 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2018.root");
	  // f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2018_Transfer_Factors.root","recreate");
	}
      else
	{ cout<<"please give a proper input"; return 0;
	}
    }

  TH1D *fail_accept,*fail_id,*fail_iso,*cr,*fake_photon,*lost_events, *srbins_lostel,*srbins_pred;
  
  fake_photon = (TH1D*)f1->Get("fake_photon_2");
  fail_accept = (TH1D*)f1->Get("fail_accept_2");
  fail_id     = (TH1D*)f1->Get("fail_id_2");
  fail_iso    = (TH1D*)f1->Get("fail_iso_2");
  cr          = (TH1D*)f1->Get("one_lep_cr_2");
  lost_events = (TH1D*)f1->Get("lost_event_2");
  // srbins_lostel = (TH1D*)f1->Get("srbins_lostel");
  // srbins_pred = (TH1D*)f1->Get("srbins_cr");
  TH1D *lost_events_copy = (TH1D*)lost_events->Clone("lost_events_copy");
  TH1D *cr_copy = (TH1D*)cr->Clone("cr_copy");
 
  //  f1->Close();
  
  
  
  // double total=0;
  //  const char* str[] = {"NJets_{0}^{=2} & 100<MET<150","NJets_{0}^{=2} & MET#geq 150","NJets_{0}^{=3} & 100<MET<150","NJets_{0}^{=3} & MET#geq 150","NJets_{0}^{=4} & 100<MET<150","NJets_{0}^{=4} & MET#geq 150","NJets_{0}^{5-6} & 100<MET<150","NJets_{0}^{5-6} & MET#geq 150","NJets_{0}^{#geq7} & 100<MET<150","NJets_{0}^{#geq7} & MET#geq 150","NJets_{= 1}^{2-4} & 100<MET<150","NJets_{= 1}^{2-4} & MET#geq 150","NJets_{= 1}^{5-6} & 100<MET<150","NJets_{= 1}^{5-6} & MET#geq 150","NJets_{= 1}^{#geq 7} & 100<MET<150","NJets_{= 1}^{#geq 7} & MET#geq 150","NJets_{#geq 2}^{2-4} & 100<MET<150","NJets_{#geq 2}^{2-4} & MET#geq 150","NJets_{#geq 2}^{5-6} & 100<MET<150","NJets_{#geq 2}^{5-6} & MET#geq 150","NJets_{#geq 2}^{#geq 7} & 100<MET<150","NJets_{#geq 2}^{#geq 7} & MET#geq 150"};

  const char* str[] = {"NJets_{0}^{=2} & 100<MET<150","NJets_{0}^{=2} & MET#geq 150","NJets_{0}^{=3} & 100<MET<150","NJets_{0}^{=3} & MET#geq 150","NJets_{0}^{=4} & 100<MET<150","NJets_{0}^{=4} & MET#geq 150","NJets_{0}^{5-6} & 100<MET<150","NJets_{0}^{5-6} & MET#geq 150","NJets_{0}^{#geq7} & 100<MET<150","NJets_{0}^{#geq7} & MET#geq 150","NJets_{#geq 1}^{2-4} & 100<MET<150","NJets_{#geq 1}^{2-4} & MET#geq 150","NJets_{#geq 1}^{5-6} & 100<MET<150","NJets_{#geq 1}^{5-6} & MET#geq 150","NJets_{#geq 1}^{#geq 7} & 100<MET<150","NJets_{#geq 1}^{#geq 7} & MET#geq 150"};


  // Making new copy of hist with total divided

  TFile *f2 = new TFile("/home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2016_Transfer_Factors.root","recreate");
  
  TH1D *total = new TH1D("total","Total = fail_id+fail_iso+fail_accept+1e_cr",16,1,17);
  for(int i=1;i<17;i++)
    { total->GetXaxis()->SetBinLabel(i,str[i-1]);}
  // total->Add(fail_accept);
  // total->Add(fail_id);
  // total->Add(fail_iso);
  total->Add(lost_events);
  total->Add(cr);
  // Making new copy of hist with total divided

  // fake_photon->Divide(total);
  // fail_accept->Divide(total);
  // fail_id->Divide(total);
  // fail_iso->Divide(total);
  lost_events->Divide(total);
  cr->Divide(total);
 
  THStack *stack = new THStack("Stack","stack hist");
  TCanvas *c1 = new TCanvas("stackhist","stackhist",1600,900);

  
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);

  TPad *pad2 = new TPad("pad1","pad1",0,0.0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad1->Draw();//pad1->SetGridx();
  pad2->Draw();pad2->SetGridx();pad2->SetGridy();
  pad1->cd();
  
  gStyle->SetPalette(kOcean);
  // fail_accept->SetFillStyle(3008);
  // fail_accept->SetFillColor(kGreen);
  // fail_id->SetFillStyle(3019);
  // fail_id->SetFillColor(kOrange);
  // fail_iso->SetFillStyle(3113);
  // fail_iso->SetFillColor(kBlue+3);
  // fake_photon->SetFillColor(kBlue-9);
  lost_events->SetFillStyle(3245);
  lost_events->SetFillColor(kViolet);
  cr->SetFillColor(kGray);
  
  TLegend *legend = new TLegend(0.6,0.8,0.9,0.9);
  stack->Add(cr);
  stack->Add(lost_events);
  // stack->Add(fail_accept);
  // stack->Add(fail_id);
  // stack->Add(fail_iso);
  //  stack->Add(fake_photon);
  stack->Draw("hist");
  legend->SetNColumns(2);
  legend->SetBorderSize(1);

  stack->GetYaxis()->SetTitleSize(0.05);
  stack->GetYaxis()->SetTitleOffset(0.85);
  stack->GetYaxis()->SetLabelSize(0.07);
  stack->GetYaxis()->SetTitle(0);
  stack->SetMaximum(1.2);
  char *title = new char[200];
  sprintf(title,"lost e^{-} -  %d",year);
  stack->SetTitle(title);

  legend->AddEntry(cr,"1L CR","f");
  legend->AddEntry(lost_events,"lost events","f");
  // legend->AddEntry(fail_accept,"Fail Accept","f");
  // legend->AddEntry(fail_id,"Fail Id","f");
  // legend->AddEntry(fail_iso,"Fail Iso","f");
  //  legend->AddEntry(fake_photon,"Fake photon","f");
  legend->SetTextSize(0.03);
  legend->Draw();


  pad2->cd();

  TH1D *TF = new TH1D("tf","Transfer factor",16,1,17);
  // for(int i=1;i<17;i++)
  //   { TF->GetXaxis()->SetBinLabel(i,str[i-1]); }
  // TF->Add(fail_accept);
  // TF->Add(fail_id);
  // TF->Add(fail_iso);
  //TF->Add(fake_photon_1);
  TF->Add(lost_events);
  TF->GetYaxis()->SetRangeUser(0,2);
  TF->Sumw2();
  TF->SetStats(0);
  TF->Divide(cr);
  TF->Draw("ep");
  TF->SetTitle(0);

  TF->GetXaxis()->SetTitle("Bin Number");
  TF->GetXaxis()->SetTitleSize(0.14);
  TF->GetXaxis()->SetTitleOffset(0.85);
  TF->GetXaxis()->SetLabelSize(0.15);

  TF->GetYaxis()->SetTitle("Transfer factor");
  TF->GetYaxis()->SetTitleOffset(0.35);
  TF->GetYaxis()->SetTitleSize(0.13);
  TF->GetYaxis()->SetLabelSize(0.09);
  TF->SetLineWidth(3);
  //TF->SetMaximum(1.99);
  //TF->SetMinimum(0.01);

  

  for(int i=1;i<=16;i++)
    { //cout<<"The Transfer Factor in bin"<<i<<" = "<<TF->GetBinContent(i)<<" +- "<<TF->GetBinError(i)<<endl;
      cout<<fixed<<setprecision(3)<<TF->GetBinContent(i)<<" +- "<<TF->GetBinError(i)<<endl;
    }

   
  // Lines and text
  pad1->cd();
  TLine *t1 = new TLine();
  t1->SetLineStyle(1);
  t1->SetLineWidth(3);
  t1->DrawLine(11,0.00,11,1.25);

  TLine *t2 = new TLine();
  t2->SetLineStyle(2);
  t2->SetLineWidth(2);
  t2->DrawLine(3,0.00,3,1.1);
  t2->DrawLine(5,0.00,5,1.1);
  t2->DrawLine(7,0.00,7,1.1);
  t2->DrawLine(9,0.00,9,1.1);
  t2->DrawLine(13,0.00,13,1.1);
  t2->DrawLine(15,0.00,15,1.1);

  TLine *t5 = new TLine();
  t5->SetLineStyle(4);
  t5->SetLineWidth(1);
  t5->DrawLine(2,0.00,2,1.0);
  t5->DrawLine(4,0.00,4,1.0);
  t5->DrawLine(6,0.00,6,1.0);
  t5->DrawLine(8,0.00,8,1.0);
  t5->DrawLine(10,0.00,10,1.0);
  t5->DrawLine(12,0.00,12,1.0);
  t5->DrawLine(14,0.00,14,1.0);
  t5->DrawLine(16,0.00,16,1.0);

  
  // Njet labels

  TLatex *text_njet = new TLatex();
  text_njet->SetTextFont(42);
  text_njet->SetTextSize(0.035);
  text_njet->SetTextAlign(22);
  text_njet->DrawLatex(2.0,1.05,"NJets^{= 0} = 2");
  text_njet->DrawLatex(4.0,1.05,"NJets^{ = 0} = 3");
  text_njet->DrawLatex(6.0,1.05,"NJets^{ = 0} = 4");
  text_njet->DrawLatex(8.0,1.05,"5 #leq NJets^{ = 0} #leq 6");
  text_njet->DrawLatex(10.0,1.05,"NJets^{ = 0} #geq 7");
  text_njet->DrawLatex(12.0,1.05,"2 #leq NJets^{ #geq 1} #leq 4");
  text_njet->DrawLatex(14.0,1.05,"5 #leq NJets^{ #geq 1} #leq 6");
  text_njet->DrawLatex(16.0,1.05,"NJets^{ #geq 1} #geq 7");

  TLatex *pt_text = new TLatex();
  pt_text->SetTextFont(42);
  pt_text->SetTextSize(0.018);
  pt_text->SetTextAlign(22);

  pt_text->DrawLatex(1.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(2.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(3.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(4.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(5.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(6.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(7.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(8.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(9.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(10.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(11.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(12.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(13.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(14.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(15.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(16.5,0.4,"p_{T}^{miss} #geq 150");
  //pt_text->DrawLatex(17.0,0.4,"100 #leq p_{T}^{miss} #leq 150");
  
  
  pad2->cd();
  TLine *t3 = new TLine();
  t3->SetLineStyle(1);
  t3->SetLineWidth(3);
  t3->DrawLine(11,0.00,11,2.0);

  TLine *t4 = new TLine();
  t4->SetLineStyle(2);
  t4->SetLineWidth(2);
  t4->DrawLine(3,0.00,3,2.0);
  t4->DrawLine(5,0.00,5,2.0);
  t4->DrawLine(7,0.00,7,2.0);
  t4->DrawLine(9,0.00,9,2.0);
  t4->DrawLine(13,0.00,13,2.0);
  t4->DrawLine(15,0.00,15,2.0);

  TLine *t6 = new TLine();
  t6->SetLineStyle(4);
  t6->SetLineWidth(1);
  t6->DrawLine(2,0.00,2,1.0);
  t6->DrawLine(4,0.00,4,1.0);
  t6->DrawLine(6,0.00,6,1.0);
  t6->DrawLine(8,0.00,8,1.0);
  t6->DrawLine(10,0.00,10,1.0);
  t6->DrawLine(12,0.00,12,1.0);
  t6->DrawLine(14,0.00,14,1.0);
  t6->DrawLine(16,0.00,16,1.0);
  
   c1->SaveAs(save);
  //  TF->Write();
  // f2->WriteTObject(TF);
  // f2->Close();
 

  ///////////////////////////           For Closure test            ///////////////////////////////

  TH1D *cr_pred;
  cr_pred      = (TH1D*)f1->Get("cr_pred_2");

  TH1D *sr_pred = new TH1D("sr_pred","SR prediction",16,1,17);
  //sr_pred->Multiply(cr_pred,TF);   // Using CR event where only NElectron==1 cut is applied 
  sr_pred->Multiply(cr_copy,TF);  // using CR where events are filtered after the gen cuts etc.

  // // TH1D *cr_pred,*sr_pred;
  
      
  TCanvas *c2 = new TCanvas("prediction","prediction",1600,900);
  c2->cd();

  TPad *pad3 = new TPad("pad3","pad3",0,0.3,1,1);
  pad3->SetBottomMargin(0);

  TPad *pad4 = new TPad("pad4","pad4",0,0.0,1,0.3);
  pad4->SetTopMargin(0);
  pad4->SetBottomMargin(0.3);
  pad3->Draw();//pad1->SetGridx();
  pad4->Draw();pad4->SetGridx();pad4->SetGridy();
  pad3->cd();
  sr_pred->SetStats(0);
  lost_events_copy->SetStats(0);
  gStyle->SetPalette(kOcean);
  lost_events_copy->SetLineStyle(3245);
  lost_events_copy->SetLineColor(kViolet);
  sr_pred->SetLineColor(kRed);
  TLegend *legend2 = new TLegend(0.6,0.8,0.9,0.9);
  sr_pred->Draw("histe");
  lost_events_copy->Draw("histe sames");
  legend2->SetNColumns(2);
  legend2->SetBorderSize(1);
  
  sr_pred->GetYaxis()->SetTitleSize(0.05);
  sr_pred->GetYaxis()->SetTitleOffset(0.85);
  sr_pred->GetYaxis()->SetLabelSize(0.07);
  sr_pred->GetYaxis()->SetTitle(0);
  sprintf(title,"lost e^{-} -  %d  (Closure Test)",year);
  sr_pred->SetTitle(title);

  legend2->AddEntry(lost_events_copy,"SR observed MC","l");
  legend2->AddEntry(sr_pred,"SR prediction using TF","l");
  legend2->SetTextSize(0.03);
  legend2->Draw();

  pad4->cd();

  TH1D *ratio = (TH1D*)sr_pred->Clone("ratio");
  
  ratio->GetYaxis()->SetRangeUser(0,2);
  ratio->Sumw2();
  ratio->SetStats(0);
  ratio->Divide(lost_events_copy);
  ratio->Draw("ep");
  ratio->SetTitle(0);

  ratio->GetXaxis()->SetTitle("Bin Number");
  ratio->GetXaxis()->SetTitleSize(0.14);
  ratio->GetXaxis()->SetTitleOffset(0.85);
  ratio->GetXaxis()->SetLabelSize(0.15);

  ratio->GetYaxis()->SetTitle("#frac{Pred}{obs}");
  ratio->GetYaxis()->SetTitleOffset(0.35);
  ratio->GetYaxis()->SetTitleSize(0.13);
  ratio->GetYaxis()->SetLabelSize(0.09);
  ratio->SetLineWidth(3);

  pad3->cd();
  t1->SetLineStyle(1);
  t1->SetLineWidth(3);
  t1->DrawLine(11,0.00,11,1.25);

  t2->SetLineStyle(2);
  t2->SetLineWidth(2);
  t2->DrawLine(3,0.00,3,1.1);
  t2->DrawLine(5,0.00,5,1.1);
  t2->DrawLine(7,0.00,7,1.1);
  t2->DrawLine(9,0.00,9,1.1);
  t2->DrawLine(13,0.00,13,1.1);
  t2->DrawLine(15,0.00,15,1.1);

  t5->SetLineStyle(4);
  t5->SetLineWidth(1);
  t5->DrawLine(2,0.00,2,1.0);
  t5->DrawLine(4,0.00,4,1.0);
  t5->DrawLine(6,0.00,6,1.0);
  t5->DrawLine(8,0.00,8,1.0);
  t5->DrawLine(10,0.00,10,1.0);
  t5->DrawLine(12,0.00,12,1.0);
  t5->DrawLine(14,0.00,14,1.0);
  t5->DrawLine(16,0.00,16,1.0);

  
  // Njet labels

  text_njet->SetTextFont(42);
  text_njet->SetTextSize(0.035);
  text_njet->SetTextAlign(22);
  text_njet->DrawLatex(2.0,1.05,"NJets^{= 0} = 2");
  text_njet->DrawLatex(4.0,1.05,"NJets^{ = 0} = 3");
  text_njet->DrawLatex(6.0,1.05,"NJets^{ = 0} = 4");
  text_njet->DrawLatex(8.0,1.05,"5 #leq NJets^{ = 0} #leq 6");
  text_njet->DrawLatex(10.0,1.05,"NJets^{ = 0} #geq 7");
  text_njet->DrawLatex(12.0,1.05,"2 #leq NJets^{ #geq 1} #leq 4");
  text_njet->DrawLatex(14.0,1.05,"5 #leq NJets^{ #geq 1} #leq 6");
  text_njet->DrawLatex(16.0,1.05,"NJets^{ #geq 1} #geq 7");

  pt_text->SetTextFont(42);
  pt_text->SetTextSize(0.018);
  pt_text->SetTextAlign(22);

  pt_text->DrawLatex(1.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(2.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(3.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(4.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(5.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(6.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(7.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(8.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(9.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(10.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(11.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(12.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(13.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(14.5,0.4,"p_{T}^{miss} #geq 150");
  pt_text->DrawLatex(15.5,0.4,"100 #leq p_{T}^{miss} #leq 150");
  pt_text->DrawLatex(16.5,0.4,"p_{T}^{miss} #geq 150");


   c2->SaveAs(save2);
}

