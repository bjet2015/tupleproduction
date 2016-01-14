/* buildtuplemc.C - makes lightweight ntuples from HiForest files
 * change samplesfolder if necessary and subfoldernames, pthats
 * change jettree if you want to use different jet algo
 * nonoverlapping = true means events from pt,hat_bin[i-1] do not overlap 
 *                  with pt,hat_bin[i]. This reduces statistics by 10-20%, 
 *                  but allows to check if large weights corrupt errors at high pt
 */


TString samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp";


vector<TString> subfoldernames={"qcd30","qcd50","qcd80","qcd120","qcd170"}; //just a typo in the foldername
vector<int> pthats = {           30,     50,     80,     120,     220};
TString outputfilename = "mcqcd.root";

//vector<TString> subfoldernames={"b30","b50","b80","b120","b170"};
//vector<int> pthats = {           30,   50,   80,   120,   170};
//TString outputfilename = "mcb.root";

TString jettree = "ak4PFJetAnalyzer/t";

//non-overlapping pt,hat bins loose 10-20% of statistics, but don't have large weights issue
bool nonoverlapping = false;


int getind (float pthat)
{
  if (pthat<pthats[0])
    cout<<"pthat "<<pthat<<" is less than minimum "<<pthats[0];

  for (int i=1;i<pthats.size();i++) 
    if (pthat<pthats[i]) return i-1;

  return pthats.size()-1;
  
}

void buildtuplemc()
{
  int bins = pthats.size();

  vector<double> weights (bins);//, x0(bins), x1(bins);
  vector<double> eventsCurbin(pthats.size());
  vector<double> eventsNextbin(pthats.size());
  vector<double> abovethresh(bins);
  weights[0] = 1;

  cout<<"Calculating weights"<<endl;
  cout<<" pthat "<<pthats[0]<<"\t"<<flush;
  for (int i=1;i<pthats.size();i++) {
    cout<<" pthat "<<pthats[i]<<"\t"<<flush;
    TFile *f0 = new TFile(Form("%s/%s/merged_HiForestAOD.root",samplesfolder.Data(),subfoldernames[i].Data()));
    TTree *t0 = (TTree *)f0->Get(jettree);
    double x0 = t0->GetEntries(Form("pthat>%d",pthats[i]));

    
    int minoverlap = nonoverlapping ? i-1 : 0; 
    double x1=0;
    for (int j=minoverlap;j<i;j++) {
      TFile *f1 = new TFile(Form("%s/%s/merged_HiForestAOD.root",samplesfolder.Data(),subfoldernames[j].Data()));
      TTree *t1 = (TTree *)f1->Get(jettree);
      x1+=t1->GetEntries(Form("pthat>%d",pthats[i]));
      f1->Close();
    }

    weights[i] = nonoverlapping ? weights[i-1]*x1/x0 : weights[i-1]*x1/(x0+x1);
    
    f0->Close();
  }

  double sumw = 0;
  for (auto w:weights) sumw+=w;
  for (int i=0;i<weights.size();i++) weights[i]/=sumw;

  cout<<endl<<"Weights : "<<endl;
  for (auto w:weights) cout<<w<<"\t";
  cout<<endl;

  int totentries = 0;

  //now fill histos
  TFile *fout = new TFile(outputfilename,"recreate");
  TNtuple *nt = new TNtuple("nt","nt","pthat:weight:dijet:genpt0:genpt1:jtpt0:jtpt1:jtphi0:jtphi1:jteta0:jteta1:discr_csvSimple0:discr_csvSimple1:refparton_flavorForB0:refparton_flavorForB1");
  
  for (int i=0;i<pthats.size();i++) {
    TString filename = TString::Format("%s/%s/merged_HiForestAOD.root",samplesfolder.Data(),subfoldernames[i].Data());
    cout<<endl<<"Processing file "<<filename<<endl;

    TFile *f = new TFile(filename);
    TTreeReader reader(jettree,f);
    TTreeReaderValue<float> pthat(reader, "pthat");
    TTreeReaderValue<int> nref(reader, "nref");
    TTreeReaderArray<float> genpt(reader, "genpt");
    TTreeReaderArray<float> jtpt(reader, "jtpt");
    TTreeReaderArray<float> jteta(reader, "jteta");
    TTreeReaderArray<float> jtphi(reader, "jtphi");
    TTreeReaderArray<float> discr_csvSimple(reader, "discr_csvSimple");
    TTreeReaderArray<int> refparton_flavorForB(reader, "refparton_flavorForB");
    
    int nev = reader.GetEntries(true);
    totentries+=nev;
    int onep = nev/100;
    int evCounter = 0;
    TTimeStamp t0;

    while (reader.Next()) {
      evCounter++;
      if (evCounter%onep==0) {
	std::cout << std::fixed;
	TTimeStamp t1; 
	cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
      }

      int ind[2];
      bool foundLJ=false, foundSJ = false;
      
      for (int j=0;j<*nref;j++) {
	if (!foundLJ) //looking for the leading jet
	  if (abs(jteta[j])<2) {
	    ind[0] = j;
	    foundLJ=true;
	    continue;
	  }
	if (foundLJ) {
	  double deltaphi = acos(cos(jtphi[j]-jtphi[ind[0]]));
	  if (abs(jteta[j])<2) { // && deltaphi>2./3.*3.142) { ----- delta-phi selection has to be made later!
	    ind[1] = j;
	    foundSJ = true;
	    break;
	  }
	}
      }

      if (nonoverlapping && getind(*pthat)!=i) continue;

      if (foundSJ && !foundLJ) cout<<"foundSJ && !foundLJ"<<endl;

      vector<float> v;

      if (foundLJ && foundSJ)
	v = {*pthat, (float)weights[getind(*pthat)], 1, //1 = dijet
	     genpt[ind[0]], genpt[ind[1]],
	     jtpt[ind[0]], jtpt[ind[1]],
	     jtphi[ind[0]], jtphi[ind[1]],
	     jteta[ind[0]], jteta[ind[1]],
	     discr_csvSimple[ind[0]], discr_csvSimple[ind[1]],
	     (float)refparton_flavorForB[ind[0]], (float)refparton_flavorForB[ind[1]] };
      else if (foundLJ && !foundSJ) 
	v = {*pthat, (float)weights[getind(*pthat)], 0, //0 = monojet
	     genpt[ind[0]], 0,
	     jtpt[ind[0]], 0,
	     jtphi[ind[0]], 0,
	     jteta[ind[0]], 0,
	     discr_csvSimple[ind[0]], 0,
	     (float)refparton_flavorForB[ind[0]], 0 };
      
      if (v.size()>0)
	nt->Fill(&v[0]);
    }
    f->Close();
  }
  
  fout->cd();
  nt->Write();
  fout->Close();

  cout<<endl;
  cout<<"Total input entries "<<totentries<<endl;

  }
