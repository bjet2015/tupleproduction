//for root headers...
//#include "parsecode.h"

TF1 *vertexWeighting(TString dtinputfile, TString mcinputfile)
{
  TFile *fmc = new TFile(mcinputfile);
  auto ntmc = (TTree *)fmc->Get("nt");
  TFile *fdt = new TFile(dtinputfile);
  auto ntdt = (TTree *)fdt->Get("nt");

  auto vdt = new TH1F("vdt","vdt",200,-15,15); vdt->Sumw2();
  auto vmc = new TH1F("vmc","vmc",200,-15,15); vmc->Sumw2();

  ntdt->Project("vdt","vz","weight*(jtpt1>100 && jtpt2>40 && dphi21>2.)");
  ntmc->Project("vmc","vz","pthatweight*(jtpt1>100 && jtpt2>40 && dphi21>2.)");
  
  vdt->Scale(1/vdt->Integral());
  vmc->Scale(1/vmc->Integral());

  TF1 *f1 = new TF1("f1", "gaus", -15, 15);
  vdt->Fit("f1", "R");

  TF1 *f2 = new TF1("f2", "gaus", -15, 15);
  vmc->Fit("f2", "R");
  
  TF1 *f = new TF1("fv","f1/f2", -15, 15);
// 
  // vdt->Rebin(4);
  // vmc->Rebin(4);

  auto vmcw = (TH1F *)vmc->Clone("vmcw");
  vmcw->Multiply(f);


  TCanvas *ccccc = new TCanvas("cc","cc",600,600);
  vdt->SetLineColor(kRed);
  vmc->SetLineColor(kBlue);
  vmcw->SetLineColor(kGreen);

  vdt->Draw();
  vmc->Draw("same");
  vmcw->Draw("same");
  //vwb->Draw();
  cout<<vdt->GetEntries()<<endl;
  cout<<vmc->GetEntries()<<endl;

  return f;

}

//there are 3 ways to do that:
//fit ratio with function
//fit data with function and return histogram of fdata/hmc
//don't fit
//(one can also fit MC, but it is crazy)
TF1 * centrWeighting(TString dtinputfile, TString mcinputfile)
{
  cout<<mcinputfile<<" "<<dtinputfile<<endl;

  TFile *fmc = new TFile(mcinputfile);
  auto ntmc = (TTree *)fmc->Get("nt");
  TFile *fdt = new TFile(dtinputfile);
  auto ntdt = (TTree *)fdt->Get("nt");
  
  int nbins=200;

  auto cdt = new TH1F("cdt","cdt",nbins,0,200); cdt->Sumw2();
  auto cmc = new TH1F("cmc","cmc",nbins,0,200); cmc->Sumw2();
  auto cwb = new TH1F("hCentrWeight","hCentrWeight",nbins,0,200); cwb->Sumw2(); cwb->SetLineColor(kRed);

  ntdt->Project("cdt","bin","weight*(jtpt1>100 && jtpt2>40 && dphi21>2. && hiHF<5500)");
  ntmc->Project("cmc","bin","pthatweight*(jtpt1>100 && jtpt2>40 && dphi21>2.)");// && pthatsample>30
  
  cdt->Scale(1/cdt->Integral());
  cmc->Scale(1/cmc->Integral());



  TCanvas *cccc = new TCanvas("c","c",600,600);
  cdt->Draw();
  cmc->Draw("same");
  cout<<cdt->GetEntries()<<endl;
  cout<<cmc->GetEntries()<<endl;


  cwb->Divide(cdt,cmc,1.,1.);

  TCanvas *ccc = new TCanvas("cc","cc",600,600);
  cwb->Draw();

  vector<float> res;
  for (int i=0;i<200;i++) res.push_back((float)(cwb->GetBinContent(i+1)));

 //in principle can use function fit
  TF1 *funcratio = new TF1("cbinratio","[0]*exp([1]*x+[2]*x*x+[3]*x*x*x)",0,200);
  // funcratio->SetParameter(0, 3.10026e+00);
  // funcratio->SetParameter(1,-4.08580e-02);
  // funcratio->SetParameter(2, 2.87782e-04);
  // funcratio->SetParameter(3,-1.54032e-06);

  funcratio->SetParameter(0, 3.10026e+00);
  funcratio->SetParameter(1,-4.08580e-02);
  funcratio->SetParameter(2, 2.87782e-04);
  funcratio->SetParameter(3,-1.54032e-06);
  // funcratio->SetParLimits(0,0,10);
  // funcratio->SetParLimits(1,-.1,0);
  // funcratio->SetParLimits(2,-.1,0);
  // funcratio->SetParLimits(3,-.1,0);


  cwb->Fit(funcratio,"Q");

//  cwb->Draw();

     return funcratio;

  
}


vector<float> centrWeightingold(TString dtinputfile, TString mcinputfile)
{
  cout<<mcinputfile<<" "<<dtinputfile<<endl;

  TFile *fmc = new TFile(mcinputfile);
  auto ntmc = (TTree *)fmc->Get("nt");
  TFile *fdt = new TFile(dtinputfile);
  auto ntdt = (TTree *)fdt->Get("nt");
  
  int nbins=200;

  auto cdt = new TH1F("cdt","cdt",nbins,0,200); cdt->Sumw2();
  auto cmc = new TH1F("cmc","cmc",nbins,0,200); cmc->Sumw2();
  auto cwb = new TH1F("hCentrWeight","hCentrWeight",nbins,0,200); cwb->Sumw2(); cwb->SetLineColor(kRed);

  ntdt->Project("cdt","bin","weight*(jtpt1>100 && jtpt2>40)");
  ntmc->Project("cmc","bin","pthatweight*(jtpt1>100 && jtpt2>40)");// && pthatsample>30
  
  cdt->Scale(1/cdt->Integral());
  cmc->Scale(1/cmc->Integral());



  TCanvas *cccc = new TCanvas("c","c",600,600);
  cdt->Draw();
  cmc->Draw("same");
  cout<<cdt->GetEntries()<<endl;
  cout<<cmc->GetEntries()<<endl;


  cwb->Divide(cdt,cmc,1.,1.);

  TCanvas *ccc = new TCanvas("cc","cc",600,600);
  cwb->Draw();

  vector<float> res;
  for (int i=0;i<200;i++) res.push_back((float)(cwb->GetBinContent(i+1)));

 //in principle can use function fit
  TF1 *funcratio = new TF1("cbinratio","[0]*exp([1]*x+[2]*x*x+[3]*x*x*x)",0,200);
  // funcratio->SetParameter(0, 3.10026e+00);
  // funcratio->SetParameter(1,-4.08580e-02);
  // funcratio->SetParameter(2, 2.87782e-04);
  // funcratio->SetParameter(3,-1.54032e-06);

  funcratio->SetParameter(0, 3.10026e+00);
  funcratio->SetParameter(1,-4.08580e-02);
  funcratio->SetParameter(2, 2.87782e-04);
  funcratio->SetParameter(3,-1.54032e-06);
  // funcratio->SetParLimits(0,0,10);
  // funcratio->SetParLimits(1,-.1,0);
  // funcratio->SetParLimits(2,-.1,0);
  // funcratio->SetParLimits(3,-.1,0);


  cwb->Fit(funcratio,"Q");

//  cwb->Draw();

     return res;

  


  // //Use the easiest way for now
  // TFile *fout = new TFile(centrfile,"recreate");
  // fout->cd();
  // cwb->Write();
  // fout->Close();


  //this is ok, but... not needed

  TCanvas *ccomp = new TCanvas("ccomp","ccomp",600,600);
  cdt->Draw();
  cmc->Draw("same");

  
  //in principle can use function fit
  TF1 *funcdt = new TF1("cbindt","[0]*exp([1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x)",0,200);
  //  funcdt->FixParameter(0,cdt->GetBinContent(1));
  funcdt->SetParameter(0,0.1);
  funcdt->SetParameter(1,-0.1);
  funcdt->SetParameter(2,1E-3);
  funcdt->SetParameter(3,-1E-5);
  funcdt->SetParameter(4,1E-8);

  cdt->Fit(funcdt,"Q");

  TCanvas *c = new TCanvas("c","c",600,600);
  cdt->Draw();
  //  fdt->Draw();

  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  cmc->Draw();
  auto cw = new TH1F("cw","cw",200,0,200); cw->Sumw2();
  for (int i=1;i<200;i++) {cw->SetBinContent(i,1); cw->SetBinError(i,0); }

  cw->Multiply(funcdt);
  cw->Divide(cmc);
  //  cmc->Divide(funcdt,cmc);
  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  cwb->Draw();
  //  cw->Draw();
  //    
 // vector<float> res;
 //   for (int i=0;i<200;i++) res.push_back((float)(cwb->GetBinContent(i+1)));

 //     return res;

  

}


//test
void weighting()
{
  //vertexWeighting("/data_CMS/cms/lisniak/bjet2015/dtPbj40akVs4PF_djt.root","/data_CMS/cms/lisniak/bjet2015/mcPbqcdakVs4PF_djt.root");
  centrWeighting("/data_CMS/cms/lisniak/bjet2015/dtPbj40akPu4PF_djt.root","/data_CMS/cms/lisniak/bjet2015/mcPbqcdakPu4PF_djt.root");
   
}
