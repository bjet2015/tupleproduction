/* buildtuplemc.C - makes lightweight ntuples from HiForest files
 * change samplesfolder if necessary and subfoldernames, pthats
 * change jettree if you want to use different jet algo
 * nonoverlapping = true means events from pt,hat_bin[i-1] do not overlap 
 *                  with pt,hat_bin[i]. This reduces statistics by 10-20%, 
 *                  but allows to check if large weights corrupt errors at high pt
 */

TString samplesfolder, outputfilename, jettree;
vector<TString> subfoldernames;
vector<int> pthats;
vector<double> CS;

void Init(TString colType, TString mcType)
{
  TString outputfolder = "/data_CMS/cms/lisniak/bjet2015/";
  if (colType=="PbPb" && mcType=="qcd") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/v2";
    subfoldernames={"qcd30","qcd50","qcd80","qcd120","qcd170"};
    pthats = {           30,     50,     80,     120,     170};
    outputfilename = outputfolder+"mcPbPbqcd_inc.root";
    jettree = "akPu4PFJetAnalyzer/t";
  }

  else  if (colType=="PbPb" && mcType=="bjet") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/v2";
    subfoldernames={"b30","b50","b80","b120","b170"};
    pthats = {           30,     50,     80,     120,     170};
    outputfilename = outputfolder+"mcPbPbbjet_inc.root";
    jettree = "akPu4PFJetAnalyzer/t";
  }

  else  if (colType=="pp" && mcType=="qcd") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/v2";
    subfoldernames={"qcd30","qcd50","qcd80","qcd120"};//,"qcd170"}; //just a typo in the foldername
    pthats = {           30,     50,     80,     120};//,     220};
    outputfilename = outputfolder+"mcppqcd_inc.root";
    jettree = "ak4PFJetAnalyzer/t";
    CS = {3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05};
  } 

  else  if (colType=="pp" && mcType=="bjet") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/v2";
    subfoldernames={"b30","b50","b80","b120","b170"};
    pthats = {           30,   50,   80,   120,   170};
    outputfilename = outputfolder+"mcppbjet_inc.root";
    jettree = "ak4PFJetAnalyzer/t";
  }
  else cout<<"Unknown type of collision-mc. Use PbPb/pp-qcd/bjet"<<endl;
  
}


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

TTree *GetTree(TFile *f, TString treename)
{
  TTree *t = (TTree *)f->Get(treename);
  //PbPb bjet pthat120 has only unsubtracted jets!!!!! 
  //TODO: figure out
  if (t==0) t = (TTree *)f->Get("ak4PFJetAnalyzer/t");
  return t;
}

vector<TString> list_files(const char *dirname="/data_CMS/cms/mnguyen/bJet2015/mc/pp/v2", const char *ext=".root")
{
  vector<TString> names;// = new vector<TString>();                                                                          
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)) {
        names.push_back(TString(dirname)+"/"+fname);
      }
    }
  }

  return names;
}


double geteventsinpthat(TString subfoldername, float pthat) {
  double x0=0;
  auto files = list_files(TString::Format("%s/%s/",samplesfolder.Data(),subfoldername.Data()));
  for (auto f:files) {
    TFile *f0 = new TFile(f);
    TTree *t0 = GetTree(f0,jettree);
    x0 += t0->GetEntries(Form("pthat>%f",pthat));
    f0->Close();
  }

  return x0;
}

double geteventsinpthat(TString subfoldername, float pthatmin, float pthatmax) {
  double x0=0;
  auto files = list_files(TString::Format("%s/%s/",samplesfolder.Data(),subfoldername.Data()));
  for (auto f:files) {
    TFile *f0 = new TFile(f);
    TTree *t0 = GetTree(f0,jettree);
    x0 += t0->GetEntries(Form("pthat>%f && pthat<%f",pthatmin, pthatmax));
    f0->Close();
  }

  return x0;
}

void buildtuplemc_inc(TString colType="PbPb", TString mcType="qcd")
{
  Init(colType, mcType);

  int bins = pthats.size();

  vector<double> weights (bins);
  vector<double> eventsCurbin(pthats.size());
  vector<double> eventsNextbin(pthats.size());
  vector<double> abovethresh(bins);

  weights[0] = 1;
  
  cout<<"Calculating weights"<<endl;
  //  cout<<" pthat "<<pthats[0]<<"\t"<<flush;
  /*  for (int i=1;i<pthats.size();i++) {
    cout<<" pthat "<<pthats[i]<<"\t"<<flush;
    double x0 = geteventsinpthat(subfoldernames[i], pthats[i]);

    
    int minoverlap = nonoverlapping ? i-1 : 0; 
    double x1=0;
    for (int j=minoverlap;j<i;j++)
      x1+=geteventsinpthat(subfoldernames[j], pthats[i]);
    
    //    cout<<pthats[i]<<" "<<x1<<" "<<x0<<" "<<weights[i-1]<<endl;
    weights[i] = nonoverlapping ? weights[i-1]*x1/x0 : weights[i-1]*x1/(x0+x1);
  }
  double sumw = 0;
  for (auto w:weights) sumw+=w;
  for (int i=0;i<weights.size();i++) weights[i]/=sumw;
  
  */

  for (int i=0;i<pthats.size();i++) {
    cout<<" pthat "<<pthats[i]<<"\t"<<flush;  
    double x1=0;
    if (i!=pthats.size()-1) {
      for (int j=0;j<=i;j++)
	x1+=geteventsinpthat(subfoldernames[j], pthats[i], pthats[i+1]);
      weights[i] = (CS[i] - CS[i+1])/x1;
    } else {
      for (int j=0;j<=i;j++)
	x1+=geteventsinpthat(subfoldernames[j], pthats[i]);
      weights[i] = CS[i]/x1;
    }
  }

  cout<<endl<<"Weights : "<<endl;
  for (auto w:weights) cout<<w<<"\t";  cout<<endl;


  int totentries = 0;

  //now fill histos
  TFile *fout = new TFile(outputfilename,"recreate");
  TNtuple *nt = new TNtuple("nt","nt","pthat:pthatbin:weight:genpt:jtpt:jtphi:jteta:discr_csvSimple:refparton_flavorForB");
  
  for (int i=0;i<pthats.size();i++) {
    auto files = list_files(TString::Format("%s/%s/",samplesfolder.Data(),subfoldernames[i].Data()));

    for (auto filename:files) {
    cout<<endl<<"Processing file "<<filename<<endl;

    TFile *f = new TFile(filename);
    //PbPb bjet pthat120 has only unsubtracted jets!!!!
    TString treename = f->Get(jettree) != 0 ? jettree : "ak4PFJetAnalyzer";
    TTreeReader reader(treename,f);
    TTreeReaderValue<float> pthat(reader, "pthat");
    TTreeReaderValue<int> nref(reader, "nref");
    TTreeReaderArray<float> genpt(reader, "genpt");
    TTreeReaderArray<float> jtpt(reader, "jtpt");
    TTreeReaderArray<float> jteta(reader, "jteta");
    TTreeReaderArray<float> jtphi(reader, "jtphi");
    TTreeReaderArray<float> discr_csvSimple(reader, "discr_csvSimple");
    TTreeReaderArray<int> refparton_flavorForB(reader, "refparton_flavorForB");

    TTreeReader readerevt("hiEvtAnalyzer/HiTree",f);
    TTreeReaderValue<float> vz(readerevt, "vz");
    
    int nev = reader.GetEntries(true);
    totentries+=nev;
    int onep = nev/100;
    int evCounter = 0;
    TTimeStamp t0;

    while (reader.Next()) {
      readerevt.Next();
      evCounter++;
      if (evCounter%onep==0) {
	std::cout << std::fixed;
	TTimeStamp t1; 
	cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
      }

      if (abs(*vz)>15) continue;

      int ind[2];
      bool foundLJ=false, foundSJ = false;
      
      for (int j=0;j<*nref;j++) {

      vector<float> v;
      v = {*pthat, (float)i, (float)weights[getind(*pthat)],
	   genpt[j], jtpt[j], jtphi[j], jteta[j], discr_csvSimple[j],
	   (float)refparton_flavorForB[j]};
      
      nt->Fill(&v[0]);
      }
    }
    f->Close();
    }
  }
  
  fout->cd();
  nt->Write();
  fout->Close();

  cout<<endl;
  cout<<"Total input entries "<<totentries<<endl;

  }
