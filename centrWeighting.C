void centrWeighting()
{
  TFile *fmc = new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbPbqcdakVs4PFJetAnalyzer_evt.root");
  auto ntmc = (TTree *)fmc->Get("nt");
  TFile *fdt = new TFile("/data_CMS/cms/lisniak/bjet2015/datamergedPbPb_akVs4PFJetAnalyzer_evt.root");
  auto ntdt = (TTree *)fdt->Get("nt");
  
  auto cdt = new TH1F("cdt","cdt",200,0,200); cdt->Sumw2();
  auto cmc = new TH1F("cmc","cmc",200,0,200); cmc->Sumw2();
  auto cwb = new TH1F("hCentrWeight","hCentrWeight",200,0,200); cwb->Sumw2(); cwb->SetLineColor(kRed);
  auto cw = new TH1F("cw","cw",200,0,200); cw->Sumw2();

  ntdt->Project("cdt","bin","weight");
  ntmc->Project("cmc","bin","pthatweight");
  
  cdt->Scale(1/cdt->Integral());
  cmc->Scale(1/cmc->Integral());

  cwb->Divide(cdt,cmc,1.,1.);

  //Use the easiest way for now
  TFile *fout = new TFile("centralityWeights.root","recreate");
  cwb->Write();
  fout->Close();

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

  cdt->Fit(funcdt);

  TCanvas *c = new TCanvas("c","c",600,600);
  cdt->Draw();
  //  fdt->Draw();

  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  cmc->Draw();

  for (int i=1;i<200;i++) {cw->SetBinContent(i,1); cw->SetBinError(i,0); }

  cw->Multiply(funcdt);
  cw->Divide(cmc);
  //  cmc->Divide(funcdt,cmc);
  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  cwb->Draw();
  //  cw->Draw();
  //    


  

}
