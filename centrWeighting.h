vector<float> centrWeighting(TString dtinputfile, TString mcinputfile)
{
  cout<<mcinputfile<<" "<<dtinputfile<<endl;

  TFile *fmc = new TFile(mcinputfile);
  auto ntmc = (TTree *)fmc->Get("nt");
  TFile *fdt = new TFile(dtinputfile);
  auto ntdt = (TTree *)fdt->Get("nt");
  
  auto cdt = new TH1F("cdt","cdt",200,0,200); cdt->Sumw2();
  auto cmc = new TH1F("cmc","cmc",200,0,200); cmc->Sumw2();
  auto cwb = new TH1F("hCentrWeight","hCentrWeight",200,0,200); cwb->Sumw2(); cwb->SetLineColor(kRed);
  auto cw = new TH1F("cw","cw",200,0,200); cw->Sumw2();

  ntdt->Project("cdt","bin","weight*(jtpt0>120 && jtpt1>30)");
  ntmc->Project("cmc","bin","pthatweight*(jtpt0>120 && jtpt1>30 && pthatbin>0)");
  
  cdt->Scale(1/cdt->Integral());
  cmc->Scale(1/cmc->Integral());

  cwb->Divide(cdt,cmc,1.,1.);
  cout<<"1"<<endl;

  TCanvas *cccc = new TCanvas("c","c",600,600);
  cdt->Draw();
  cmc->Draw("same");
  cout<<cdt->GetEntries()<<endl;
  cout<<cmc->GetEntries()<<endl;


 vector<float> res;
   for (int i=0;i<200;i++) res.push_back((float)(cwb->GetBinContent(i+1)));

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
 // vector<float> res;
 //   for (int i=0;i<200;i++) res.push_back((float)(cwb->GetBinContent(i+1)));

 //     return res;

  

}
