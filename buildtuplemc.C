/* buildtuplemc.C - makes lightweight ntuples from HiForest files
 * change samplesfolder if necessary and subfoldernames, pthats
 * change jettree if you want to use different jet algo
 * nonoverlapping = true means events from pt,hat_bin[i-1] do not overlap 
 *                  with pt,hat_bin[i]. This reduces statistics by 10-20%, 
 *                  but allows to check if large weights corrupt errors at high pt
 */

#include <fstream>
#include "parsecode.h"

TString samplesfolder, jettree;
vector<TString> subfoldernames;
vector<int> pthats;
vector<double> CS;

TString outputfolder = "/data_CMS/cms/lisniak/bjet2015/";

void Init(bool PbPb, TString sample)
{
  if (PbPb && sample=="qcd") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"qcd30","qcd50","qcd80","qcd120"};
    pthats = {           30,     50,     80,     120};
  }
  else if (PbPb && sample=="qp8") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia8";
    subfoldernames={"qcd30","qcd50","qcd80","qcd120"};
    pthats = {           30,     50,     80,     120};
  }
  else  if (PbPb && sample=="bjt") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"b30","b50","b80","b120","b170"};
    pthats = {           30,     50,     80,     120,     170};
  }
  else  if (!PbPb && sample=="qp8") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia8";
    subfoldernames={"qcd30","qcd50","qcd80","qcd120"};
    pthats = {           30,     50,     80,     120};
    CS = {3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05};
  } 
  else  if (!PbPb && sample=="qcd") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia6";
    subfoldernames={"qcd30","qcd50","qcd80","qcd120"};
    pthats = {           30,     50,     80,     120};
    CS = {3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05};
  } 
  else  if (!PbPb && sample=="bjt") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia8";
    subfoldernames={"b30","b50","b80","b120","b170"};
    pthats = {           30,   50,   80,   120,   170};
  }
  else cout<<"Unknown type of collision-mc. Use PbPb/pp-qcd/bjet"<<endl;
}

TString getfoldername(int binnumber)
{
  return TString::Format("%s/%s/merged_HiForestAOD.root",samplesfolder.Data(),subfoldernames[binnumber].Data());
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
  //  if (t==0) {
    //    t = (TTree *)f->Get("ak3PFJetAnalyzer/t");
  //    cout<<"tree not found, using ak3PFJetAnalyzer"<<endl;
  //  }
  return t;
}

vector<TString> list_files(const char *dirname, const char *ext=".root")
{
  vector<TString> names;
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
  if (names.size()==0) return {dirname};

  return names;
}
double geteventsinpthat(int binnumber, float pthatmin)//TString subfoldername, float pthat) {
{
  double x0=0;
  auto files = list_files(getfoldername(binnumber));
  for (auto f:files) {
    TFile *f0 = TFile::Open(f);
    TTree *t0 = GetTree(f0,jettree);
    x0 += t0->GetEntries(Form("pthat>%f",pthatmin));
    f0->Close();
  }

  return x0;
}

double geteventsinpthat(int binnumber, float pthatmin, float pthatmax)
{
  double x0=0;
  auto files = list_files(getfoldername(binnumber));
  for (auto f:files) {
    TFile *f0 = TFile::Open(f);
    TTree *t0 = GetTree(f0,jettree);
    x0 += t0->GetEntries(Form("pthat>%f && pthat<%f",pthatmin, pthatmax));
    f0->Close();
  }
  return x0;
}

vector<double> calculateWeights()
{
  int bins = pthats.size();
  vector<double> weights (bins);
  weights[0] = 1;

  cout<<"Calculating weights"<<endl;
  //with distributions matching
  cout<<" pthat "<<pthats[0]<<"\t"<<flush;
  for (int i=1;i<pthats.size();i++) {
    cout<<" pthat "<<pthats[i]<<"\t"<<flush;

    TFile *f0 = new TFile(Form("%s/%s/merged_HiForestAOD.root",samplesfolder.Data(),subfoldernames[i].Data()));
    TTree *t0 = GetTree(f0,jettree);
    double x0 = t0->GetEntries(Form("pthat>%d",pthats[i]));

    int minoverlap = nonoverlapping ? i-1 : 0;
    double x1=0;
    for (int j=minoverlap;j<i;j++) {
      TFile *f1 = new TFile(Form("%s/%s/merged_HiForestAOD.root",samplesfolder.Data(),subfoldernames[j].Data()));
      TTree *t1 = GetTree(f1,jettree);
      x1+=t1->GetEntries(Form("pthat>%d",pthats[i]));
      f1->Close();
    }

    weights[i] = nonoverlapping ? weights[i-1]*x1/x0 : weights[i-1]*x1/(x0+x1);

    f0->Close();
  }
  /*
  //with cross-sections
  for (int i=0;i<pthats.size();i++) {
    cout<<" pthat "<<pthats[i]<<"\t"<<flush;
    double x1=0;
    if (i!=pthats.size()-1) {
      for (int j=0;j<=i;j++)
        x1+=geteventsinpthat(i, pthats[i], pthats[i+1]);
      weights[i] = (CS[i] - CS[i+1])/x1;
    } else {
      for (int j=0;j<=i;j++)
        x1+=geteventsinpthat(i, pthats[i]);
      weights[i] = CS[i]/x1;
    }
  }
  */
  

  double sumw = 0;
  for (auto w:weights) sumw+=w;
  for (int i=0;i<weights.size();i++) weights[i]/=sumw;

  cout<<endl<<"Weights : "<<endl;
  for (auto w:weights) cout<<w<<"\t";
  cout<<endl;


  return weights;
}

bool B(int j) { return abs(j)==5; }
bool C(int j) { return abs(j)==4; }
bool G(int j) { return abs(j)==21; }
bool L(int j) { return abs(j)==1 || abs(j)==2 || abs(j)==3; }
bool X(int j) { return !(B(j) || C(j));}

int getPairCode(int lj, int sj)
{
  if (B(lj) && B(sj)) return 0;
  if (C(lj) && C(sj)) return 1;
  if ((B(lj) && C(sj)) || (C(lj) && B(sj))) return 2;
  if ((B(lj) && X(sj)) || (X(lj) && B(sj))) return 3;
  if ((X(lj) && X(sj)) || (C(lj) && X(sj)) || (X(lj) && C(sj))) return 4;

  return 5;
}

bool file_exist(const char *fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}

vector<float> getCentralityWeights()
{
  TString filename = "centralityWeights.root";
  Long_t *id,*size,*flags,*mt;
  TH1F *h = 0;
  if (file_exist(filename)) {
    TFile *f = new TFile(filename);
    h = (TH1F *)f->Get("hCentrWeight");
  }
  vector<float> res;
  for (int i=0;i<200;i++) res.push_back(h!=0 ? (float)(h->GetBinContent(i)) : 0);
  return res;
}

void buildtuplemc(TString code)
{
  if (!mc(code)) { cout<<"Not mc: "<<code<<", exiting..."<<endl; return;}

  bool PbPb = isPbPb(code);
  TString sample = getSample(code);
  jettree = getjettree(code);

  Init(PbPb, sample);

  TString outputfilenamedj = outputfolder+"/"+code+"_djt.root";
  TString outputfilenameinc = outputfolder+"/"+code+"_inc.root";
  TString outputfilenameevt = outputfolder+"/"+code+"_evt.root";


  auto weights = calculateWeights();
  auto centrWeights = getCentralityWeights();
  int totentries = 0;

  //now fill histos
  TFile *foutdj = new TFile(outputfilenamedj,"recreate");
  TNtuple *ntdj = new TNtuple("nt","nt","pthat:pthatbin:pthatweight:bin:centrWeight:weight:dijet:subid0:subid1:genpt0:genpt1:rawpt0:rawpt1:jtpt0:jtpt1:jtphi0:jtphi1:jteta0:jteta1:discr_csvSimple0:discr_csvSimple1:refparton_flavorForB0:refparton_flavorForB1:pairCode");
  TFile *foutinc = new TFile(outputfilenameinc,"recreate");
  TNtuple *ntinc = new TNtuple("nt","nt","pthat:pthatbin:pthatweight:bin:centrWeight:weight:subid:genpt:rawpt:jtpt:jtphi:jteta:discr_csvSimple:refparton_flavorForB");

  TFile *foutevt = new TFile(outputfilenameevt,"recreate");
  TNtuple *ntevt = new TNtuple("nt","nt","pthat:pthatbin:pthatweight:bin:centrWeight:weight");

  for (int i=0;i<pthats.size();i++) {
    TString filename = getfoldername(i);
    cout<<endl<<"Processing file "<<filename<<endl;

    TFile *f = TFile::Open(filename);//new TFile(filename);
    //PbPb bjet pthat120 has only unsubtracted jets!!!!
    TString treename = jettree;//f->Get(jettree) != 0 ? jettree : "ak3PFJetAnalyzer";
    //    if (treename!=jettree) cout <<"Changed tree to "<<treename<<endl;
    TTreeReader reader(treename,f);
    TTreeReaderValue<float> pthat(reader, "pthat");
    TTreeReaderValue<int> nref(reader, "nref");
    TTreeReaderArray<float> genpt(reader, "genpt");
    TTreeReaderArray<float> rawpt(reader, "rawpt");
    TTreeReaderArray<float> jtpt(reader, "jtpt");
    TTreeReaderArray<float> jteta(reader, "jteta");
    TTreeReaderArray<float> jtphi(reader, "jtphi");
    TTreeReaderArray<int> subid(reader, "subid");
    TTreeReaderArray<float> discr_csvSimple(reader, "discr_csvSimple");
    TTreeReaderArray<int> refparton_flavorForB(reader, "refparton_flavorForB");
    
    TTreeReader readerevt("hiEvtAnalyzer/HiTree",f);
    TTreeReaderValue<float> vz(readerevt, "vz");
    TTreeReaderValue<int> bin(readerevt, "hiBin");

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

      int b = *bin;
      float centrWeight = PbPb ? centrWeights[b] : 1;

      vector<float> vevt = {*pthat, (float)i, (float)weights[getind(*pthat)], (float)*bin, centrWeight,
			    (float)weights[getind(*pthat)]*centrWeight};
      ntevt->Fill(&vevt[0]);


      if (abs(*vz)>15) continue;

      int ind0, ind1; //indices of leading/subleading jets in jet array
      bool foundLJ=false, foundSJ = false; //found/not found yet, for convenience

      
      for (int j=0;j<*nref;j++) {
        if (abs(jteta[j])>2) continue;

	//for inclusive plots, subid==0 everywhere
	if (subid[j]==0) {
	  vector<float> vinc;
	  vinc = {*pthat, (float)i, (float)weights[getind(*pthat)], (float)*bin, centrWeight,
		  (float)weights[getind(*pthat)]*centrWeight,
		  (float)subid[j], genpt[j], rawpt[j],jtpt[j], jtphi[j], jteta[j], discr_csvSimple[j],
		  (float)refparton_flavorForB[j]};
	  ntinc->Fill(&vinc[0]);
	}

        if (!foundLJ && subid[j]==0) { //looking for the leading jet from signal
            ind0 = j;
            foundLJ=true;
            continue;
          }
        if (foundLJ && !foundSJ) {
            ind1 = j;
            foundSJ = true;
        }
      }


      if (nonoverlapping && getind(*pthat)!=i) continue;

      vector<float> vdj;

      if (foundLJ && foundSJ)
        vdj = {*pthat,(float)i, (float)weights[getind(*pthat)], (float)*bin, centrWeight, 
	       (float)weights[getind(*pthat)]*centrWeight, 1, //1 = dijet
	       (float)subid[ind0], (float)subid[ind1],
	       genpt[ind0], genpt[ind1],
	       rawpt[ind0], rawpt[ind1],
	       jtpt[ind0], jtpt[ind1],
	       jtphi[ind0], jtphi[ind1],
	       jteta[ind0], jteta[ind1],
	       discr_csvSimple[ind0], discr_csvSimple[ind1],
	       (float)refparton_flavorForB[ind0], (float)refparton_flavorForB[ind1],
	       (float)getPairCode(refparton_flavorForB[ind0],refparton_flavorForB[ind1])};
      else if (foundLJ && !foundSJ) 
        vdj = {*pthat, (float)i, (float)weights[getind(*pthat)], (float)*bin, centrWeight,
	       (float)weights[getind(*pthat)]*centrWeight, 0, //0 = monojet
	       (float)subid[ind0], 0,
	       genpt[ind0], 0,
	       rawpt[ind0], 0,
	       jtpt[ind0], 0,
	       jtphi[ind0], 0,
	       jteta[ind0], 0,
	       discr_csvSimple[ind0], 0,
	       (float)refparton_flavorForB[ind0], 0,
	       0};
      
      if (vdj.size()>0)
        ntdj->Fill(&vdj[0]);
    }
    f->Close();
  }
  
  foutevt->cd();
  ntevt->Write();
  foutevt->Close();

  foutdj->cd();
  ntdj->Write();
  foutdj->Close();

  foutinc->cd();
  ntinc->Write();
  foutinc->Close();

  cout<<endl;
  cout<<"Total input entries "<<totentries<<endl;

  //making centrality-dependent ntuples
  PutInCbins(outputfolder, code, {{0,40}, {80,200}});

}
