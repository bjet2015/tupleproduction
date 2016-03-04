/* buildtuplemc.C - makes lightweight ntuples from HiForest files
 * change samplesfolder if necessary and subfoldernames, pthats
 * change jettree if you want to use different jet algo
 * nonoverlapping = true means events from pt,hat_bin[i-1] do not overlap 
 *                  with pt,hat_bin[i]. This reduces statistics by 10-20%, 
 *                  but allows to check if large weights corrupt errors at high pt
 */

#include <fstream>
#include "parsecode.h"
#include "centrWeighting.h"

bool PbPb = false;
bool newFlavorProcess = false;
TString samplesfolder, jettree,sample;
vector<TString> subfoldernames;
vector<int> pthats;
vector<double> CS;
vector<double> filterefficiency;
TString outputfilenamedj; 
TString outputfilenameinc;
TString outputfilenameevt;

TString outputfolder = "/data_CMS/cms/lisniak/bjet2015/";

const int NaN = -999;


void Init()
{
  filterefficiency = {1.,1.,1.,1.,1.};

  if (PbPb && sample=="qcd") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"qcd30","qcd50","qcd80","qcd120"};
    pthats = {           30,     50,     80,     120};
    CS     = {3.360e-02,3.768e-03,4.418e-04,6.159e-05};
  }
  else if (PbPb && sample=="qp8") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia8";
    subfoldernames={"qcd30","qcd50","qcd80","qcd120"};
    pthats        ={30,50,80,120};
    CS = {           3.455E-02,     4.068E-03,     4.959E-04,     7.096E-05};
  }
  else  if (PbPb && sample=="bjt") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"bjet30","bjet50","bjet80","bjet120"};
    pthats = {           30,     50,     80,     120};
    CS     = {3.360e-02,3.768e-03,4.418e-04,6.159e-05};
    filterefficiency = {6.269e-02,7.557e-02,8.640e-02,9.265e-02};
  }
  else  if (!PbPb && sample=="qp8") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia8";
    subfoldernames={"qcd30","qcd50","qcd80","qcd120"};
    pthats = {           30,     50,     80,     120};
    CS = {3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05};
  } 
  else  if (!PbPb && sample=="bp8") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia8";
    subfoldernames={"b30","b50","b80","b120"};
    pthats = {           30,     50,     80,     120};  
    CS = {3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05}; 
    filterefficiency = {8.202e-02,6.674e-02,9.647e-02,1.051e-01}; 

  } 
  else  if (!PbPb && sample=="qcd") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia6";
    subfoldernames={"qcd30","qcd50","qcd80","qcd100","qcd120"}; //"qcd100"
    pthats = {           30,     50,     80,  100,  120}; //100
    CS     = {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05}; //1.511e-04,
  } 
  else  if (!PbPb && sample=="bjt") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia6";
    subfoldernames={"bjet30","bjet50","bjet80","bjet100","bjet120"};//,"b170"};
    pthats = {           30,   50,   80,   100,     120};//,   170};
    CS     = {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
    filterefficiency = {6.269e-02,7.557e-02,8.640e-02,8.908e-02,9.265e-02};
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

vector<TString> list_files(const char *dirname, const char *exp=".*HiForestAOD.*\\.root")
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
      if (!file->IsDirectory() && fname.Contains(TRegexp(exp))) {
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
  /*
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
    }*/
  
  //with cross-sections
  /*  for (int i=0;i<pthats.size();i++) {
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
  //simpler and faster - merge spectrum, and cut onto pieces
  TH1F *hpthat = new TH1F("hpthat","hpthat",400,0,400);
  
  for (int i=0;i<pthats.size();i++) {                                                                                    
    cout<<" pthat "<<pthats[i]<<"\t"<<flush;

    auto files = list_files(TString::Format("%s/%s/",samplesfolder.Data(),subfoldernames[i].Data()));

    for (auto filename:files) {
      TFile *f = new TFile(filename);
      TH1F *hpthattemp = new TH1F("hpthattemp","hpthattemp",400,0,400);
      
      auto t = (TTree *)f->Get(jettree);

      cout<<filename<<" " <<t->GetEntries()<<endl;

      t->Project("hpthattemp","pthat");
      //      cout<<"int "<<filename<<" : "<<hpthattemp->Integral()<<endl;
      hpthat->Add(hpthattemp);
    }
  }

  for (int i=0;i<pthats.size();i++) {
    float csi = CS[i]*filterefficiency[i];

    if (i<pthats.size()-1) {
      float csi_1 = CS[i+1]*filterefficiency[i+1];
      weights[i] = (csi-csi_1)/hpthat->Integral(hpthat->FindBin(pthats[i]),hpthat->FindBin(pthats[i+1])-1);
      //      cout<<"CS : "<<csi-csi_1<<" : "<<hpthat->Integral(hpthat->FindBin(pthats[i]),hpthat->FindBin(pthats[i+1]))<<endl;
    } else {
      weights[i] = csi/hpthat->Integral(hpthat->FindBin(pthats[i]),hpthat->GetNbinsX());
      cout<<"CS : "<<csi<<" : "<<hpthat->Integral(hpthat->FindBin(pthats[i]),hpthat->GetNbinsX())<<endl;
    }
  }

  //  double sumw = 0;
  //  for (auto w:weights) sumw+=w;
  //  for (int i=0;i<weights.size();i++) weights[i]/=sumw;

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


float getFlavorProcess(TTreeReaderArray<bool> *&refparton_isGSP, TTreeReaderArray<int> *&refparton_flavorProcess, int index)
{
    //these two samples have no info about primary
  if (sample=="qp8" ||sample=="bp8") return -1;


  if (newFlavorProcess) {
    return (*refparton_flavorProcess)[index];
  } else {
    bool gsp = (*refparton_isGSP)[index];
    if (!gsp) return 1;
    else return 6;
  }
}



void do_buildtuplemc(TString code)
{
  auto weights = calculateWeights();
  //auto centrWeights = getCentralityWeights(centralityfile);
  int totentries = 0;

  //put only pthat weight
  TString djvars = TString("run:lumi:event:pthat:pthatsample:sampleEventNumber:pthatweight:bin:vz:bProdCode:dijet:bkgLeadingJet:")+
      "subid1:refpt1:rawpt1:jtpt1:jtphi1:jteta1:discr_csvSimple1:refparton_flavorForB1:refparton_flavorProcess1:svtxm1:discr_prob1:svtxdls1:svtxpt1:svtxntrk1:nsvtx1:nselIPtrk1:"+
      "subid2:refpt2:rawpt2:jtpt2:jtphi2:jteta2:discr_csvSimple2:refparton_flavorForB2:refparton_flavorProcess2:svtxm2:discr_prob2:svtxdls2:svtxpt2:svtxntrk2:nsvtx2:nselIPtrk2:dphi21:pairCode21:"+
      "subid3:refpt3:rawpt3:jtpt3:jtphi3:jteta3:discr_csvSimple3:refparton_flavorForB3:refparton_flavorProcess3:svtxm3:discr_prob3:svtxdls3:svtxpt3:svtxntrk3:nsvtx3:nselIPtrk3:dphi31:dphi32:pairCode31:pairCode32:"+
      "SLord:subidSL:refptSL:rawptSL:jtptSL:jtphiSL:jtetaSL:discr_csvSimpleSL:refparton_flavorForBSL:refparton_flavorProcessSL:svtxmSL:discr_probSL:svtxdlsSL:svtxptSL:svtxntrkSL:nsvtxSL:nselIPtrkSL:dphiSL1:pairCodeSL1";


  TString incvars = TString("run:lumi:event:pthat:pthatsample:pthatweight:bin:subid:refpt:rawpt:jtpt:jtphi:jteta:discr_csvSimple:")
      +"refparton_flavorForB:isPrimary:refparton_flavorProcess:svtxm:discr_prob:svtxdls:svtxpt:svtxntrk:nsvtx:nselIPtrk";

  //now fill histos
  TFile *foutdj = new TFile(outputfilenamedj,"recreate");
  TNtuple *ntdj = new TNtuple("nt","nt",djvars);
  TFile *foutinc = new TFile(outputfilenameinc,"recreate");
  TNtuple *ntinc = new TNtuple("nt","nt",incvars);

  TFile *foutevt = new TFile(outputfilenameevt,"recreate");
  TNtuple *ntevt = new TNtuple("nt","nt","pthat:pthatbin:pthatweight:bin");

  for (int i=0;i<pthats.size();i++) {

    auto files = list_files(TString::Format("%s/%s/",samplesfolder.Data(),subfoldernames[i].Data()));
    for (auto filename:files) {
      cout<<endl<<"Processing file "<<filename<<endl;

      //    TString filename = getfoldername(i);
      //    cout<<endl<<"Processing file "<<filename<<endl;

    TFile *f = TFile::Open(filename);//new TFile(filename);
    //PbPb bjet pthat120 has only unsubtracted jets!!!!
    TString treename = jettree;//f->Get(jettree) != 0 ? jettree : "ak3PFJetAnalyzer";
    //    if (treename!=jettree) cout <<"Changed tree to "<<treename<<endl;
    TTreeReader reader(treename,f);
    TTreeReaderValue<float> pthat(reader, "pthat");
    TTreeReaderValue<int> nref(reader, "nref");
    TTreeReaderArray<float> refpt(reader, "refpt");
    TTreeReaderArray<float> rawpt(reader, "rawpt");
    TTreeReaderArray<float> jtpt(reader, "jtpt");
    TTreeReaderArray<float> jteta(reader, "jteta");
    TTreeReaderArray<float> jtphi(reader, "jtphi");
    TTreeReaderArray<int> subid(reader, "subid");
    TTreeReaderArray<float> discr_csvSimple(reader, "discr_csvSimple");
    TTreeReaderArray<int> refparton_flavorForB(reader, "refparton_flavorForB");
    TTreeReaderArray<bool> *refparton_isGSP;

    newFlavorProcess = !PbPb && sample=="bjt";

    if (sample!="qp8" && sample!="bp8" && !newFlavorProcess)
      refparton_isGSP = new TTreeReaderArray<bool>(reader, "refparton_isGSP");
    TTreeReaderArray<int> *refparton_flavorProcess;
    TTreeReaderValue<int> *bProdCode;



    if (newFlavorProcess) {
      bProdCode = new TTreeReaderValue<int>(reader, "bProdCode");
      refparton_flavorProcess = new TTreeReaderArray<int>(reader, "refparton_flavorProcess");
    }



    TTreeReaderArray<float> discr_prob(reader, "discr_prob");
    TTreeReaderArray<float> svtxm(reader, "svtxm");
    TTreeReaderArray<float> svtxdls(reader, "svtxdls");
    TTreeReaderArray<float> svtxpt(reader, "svtxpt");
    
    TTreeReaderArray<int> svtxntrk(reader, "svtxntrk");
    TTreeReaderArray<int> nsvtx(reader, "nsvtx");
    TTreeReaderArray<int> nselIPtrk(reader, "nselIPtrk");

    TTreeReader readerevt("hiEvtAnalyzer/HiTree",f);
    TTreeReaderValue<float> vz(readerevt, "vz");
    TTreeReaderValue<int> bin(readerevt, "hiBin");
    TTreeReaderValue<unsigned int> run(readerevt, "run");
    TTreeReaderValue<unsigned int> lumi(readerevt, "lumi");
    TTreeReaderValue<unsigned long long> event(readerevt, "evt");


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
      float centrWeight = 1;//PbPb ? centrWeights[b] : 1;

      vector<float> vevt = {*pthat, (float)i, (float)weights[getind(*pthat)], (float)*bin, centrWeight,
          (float)weights[getind(*pthat)]*centrWeight};
      ntevt->Fill(&vevt[0]);


      int ind1=-1, ind2=-1, ind3=-1, indSL=-1; //indices of leading/subleading jets in jet array
      int SLord = 0;
      bool foundJ1=false, foundJ2 = false, foundJ3 = false, foundSL = false; //found/not found yet, for convenience
      bool bkgJ1 = false;

      if (abs(*vz)<15) //event-level cuts, if not passed all foundXY = false
        for (int j=0;j<*nref;j++) {
          if (abs(jteta[j])>2) continue;

          //for inclusive plots, subid==0 everywhere
          if (subid[j]==0) {
            vector<float> vinc = {(float)*run, (float)*lumi, (float)*event, *pthat, (float)pthats[i],(float)weights[getind(*pthat)],(float)*bin,
              (float)subid[j], refpt[j], rawpt[j],jtpt[j], jtphi[j], jteta[j], discr_csvSimple[j],
              (float)refparton_flavorForB[j], getFlavorProcess(refparton_isGSP,refparton_flavorProcess,j),
              svtxm[j],discr_prob[j],svtxdls[j],svtxpt[j],(float)svtxntrk[j],(float)nsvtx[j],(float)nselIPtrk[j]};

              ntinc->Fill(&vinc[0]);
            }

            //if background jumped above signal (or highest signal is outside acceptance) - then it's a "bad" event
            //this condition is propagated on every clause b/c we still need to loop jets for inclusive ntuple above
            if (!foundJ1 && subid[j]!=0) bkgJ1 = true;

          //looking for the leading jet from signal
            if (!bkgJ1 && !foundJ1 && subid[j]==0) { 
              ind1 = j;
              foundJ1=true;
            } else
            if (!bkgJ1 && foundJ1 && !foundJ2) {
              ind2 = j;
              foundJ2 = true;
            } else
            if (!bkgJ1 && foundJ1 && foundJ2 && !foundJ3) {
              ind3 = j;
              foundJ3 = true;
            }

          //we need ordinal number of SL jets, so counting until found
          //indSL != SLord because some jets are not in acceptance region
            if (!bkgJ1 && !foundSL) SLord++;

          //ind1!=j otherwise SL will be = J1
            if (!bkgJ1 && foundJ1 && ind1!=j && !foundSL && discr_csvSimple[j]>0.9) {
              indSL = j;
              foundSL = true;
            }
          }



      //if (nonoverlapping && getind(*pthat)!=i) continue;

      vector<float> vdj;

      vdj = {(float)*run, (float)*lumi, (float)*event, 
        *pthat,(float)pthats[i], (float)evCounter-1, (float)weights[getind(*pthat)], (float)*bin, *vz,
        newFlavorProcess ? (float)*(*bProdCode) : NaN,
        foundJ1 && foundJ2 ? (float)1 : (float)0, (float)bkgJ1,

        foundJ1 ? (float)subid[ind1] : NaN,
        foundJ1 ? refpt[ind1] : NaN,
        foundJ1 ? rawpt[ind1] : NaN,
        foundJ1 ? jtpt[ind1] : NaN,
        foundJ1 ? jtphi[ind1] : NaN,
        foundJ1 ? jteta[ind1] : NaN,
        foundJ1 ? discr_csvSimple[ind1] : NaN,
        foundJ1 ? (float)refparton_flavorForB[ind1] : NaN,
        foundJ1 ? getFlavorProcess(refparton_isGSP,refparton_flavorProcess,ind1) : NaN,
        foundJ1 ? svtxm[ind1] : NaN,
        foundJ1 ? discr_prob[ind1] : NaN,
        foundJ1 ? svtxdls[ind1] : NaN,
        foundJ1 ? svtxpt[ind1] : NaN,
        foundJ1 ? (float)svtxntrk[ind1] : NaN,
        foundJ1 ? (float)nsvtx[ind1] : NaN,
        foundJ1 ? (float)nselIPtrk[ind1] : NaN,
    
        foundJ2 ? (float)subid[ind2] : NaN,
        foundJ2 ? refpt[ind2] : NaN,
        foundJ2 ? rawpt[ind2] : NaN,
        foundJ2 ? jtpt[ind2] : NaN,
        foundJ2 ? jtphi[ind2] : NaN,
        foundJ2 ? jteta[ind2] : NaN,
        foundJ2 ? discr_csvSimple[ind2] : NaN,
        foundJ2 ? (float)refparton_flavorForB[ind2] : NaN,
        foundJ2 ? getFlavorProcess(refparton_isGSP,refparton_flavorProcess,ind2) : NaN,
        foundJ2 ? svtxm[ind2] : NaN,
        foundJ2 ? discr_prob[ind2] : NaN,
        foundJ2 ? svtxdls[ind2] : NaN, 
        foundJ2 ? svtxpt[ind2] : NaN,
        foundJ2 ? (float)svtxntrk[ind2] : NaN,
        foundJ2 ? (float)nsvtx[ind2] : NaN,
        foundJ2 ? (float)nselIPtrk[ind2] : NaN,
        foundJ2 && foundJ1 ? acos(cos(jtphi[ind2]-jtphi[ind1])) : NaN,
        foundJ2 && foundJ1 ? (float)getPairCode(refparton_flavorForB[ind2],refparton_flavorForB[ind1]) : NaN,
    
        foundJ3 ? (float)subid[ind3] : NaN,
        foundJ3 ? refpt[ind3] : NaN,
        foundJ3 ? rawpt[ind3] : NaN,
        foundJ3 ? jtpt[ind3] : NaN,
        foundJ3 ? jtphi[ind3] : NaN,
        foundJ3 ? jteta[ind3] : NaN,
        foundJ3 ? discr_csvSimple[ind3] : NaN,
        foundJ3 ? (float)refparton_flavorForB[ind3] : NaN,
        foundJ3 ? getFlavorProcess(refparton_isGSP,refparton_flavorProcess,ind3) : NaN,
        foundJ3 ? svtxm[ind3] : NaN,
        foundJ3 ? discr_prob[ind3] : NaN,
        foundJ3 ? svtxdls[ind3] : NaN, 
        foundJ3 ? svtxpt[ind3] : NaN,
        foundJ3 ? (float)svtxntrk[ind3] : NaN,
        foundJ3 ? (float)nsvtx[ind3] : NaN,
        foundJ3 ? (float)nselIPtrk[ind3] : NaN,
        foundJ3 && foundJ1 ? acos(cos(jtphi[ind3]-jtphi[ind1])) : NaN,
        foundJ3 && foundJ2 ? acos(cos(jtphi[ind3]-jtphi[ind2])) : NaN,
        foundJ3 && foundJ1 ? (float)getPairCode(refparton_flavorForB[ind3],refparton_flavorForB[ind1]) : NaN,
        foundJ3 && foundJ2 ? (float)getPairCode(refparton_flavorForB[ind3],refparton_flavorForB[ind2]) : NaN,

        foundSL ? (float)SLord : NaN,
        foundSL ? (float)subid[indSL] : NaN,
        foundSL ? refpt[indSL] : NaN,
        foundSL ? rawpt[indSL] : NaN,
        foundSL ? jtpt[indSL] : NaN,
        foundSL ? jtphi[indSL] : NaN,
        foundSL ? jteta[indSL] : NaN,
        foundSL ? discr_csvSimple[indSL] : NaN,
        foundSL ? (float)refparton_flavorForB[indSL] : NaN,
        foundSL ? getFlavorProcess(refparton_isGSP,refparton_flavorProcess,indSL) : NaN,
        foundSL ? svtxm[indSL] : NaN,
        foundSL ? discr_prob[indSL] : NaN,
        foundSL ? svtxdls[indSL] : NaN, 
        foundSL ? svtxpt[indSL] : NaN,
        foundSL ? (float)svtxntrk[indSL] : NaN,
        foundSL ? (float)nsvtx[indSL] : NaN,
        foundSL ? (float)nselIPtrk[indSL] : NaN,
        foundSL && foundJ1 ? acos(cos(jtphi[indSL]-jtphi[ind1])) : NaN,
        foundSL && foundJ1 ? (float)getPairCode(refparton_flavorForB[indSL],refparton_flavorForB[ind1]) : NaN

      };
      
      ntdj->Fill(&vdj[0]);
    }
    
    f->Close();
    }
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
  //  PutInCbins(outputfolder, code, {{0,40}, {80,200}});

}

void update(TString filename, vector<float> cweights)
{
  cout<<"Updating reweighting..."<<endl;
  auto f = new TFile(filename,"update");

  auto nt = (TTree *)f->Get("nt");

  float weight,cweight,pthatweight;
  float bin;
  TBranch *cw, *w;
  if (PbPb) {
    cw =  nt->Branch("centrWeight",&cweight);
    nt->SetBranchAddress("bin",&bin);
  }
  w =  nt->Branch("weight",&weight);
  nt->SetBranchAddress("pthatweight",&pthatweight);

  int n = nt->GetEntries();
  int onep = n/100;
  for (int i=0;i<n;i++) {
    //if (i%onep==0) cout<<i/onep<<endl;
    nt->GetEntry(i);


    cweight=PbPb ? cweights[bin-1] : 1;
    weight=pthatweight*cweight;

    if (PbPb)
      cw->Fill();
    w->Fill();
  }

  nt->Write();

  f->Close();


}


void buildtuplemc(TString code)
{
  if (!mc(code)) { cout<<"Not mc: "<<code<<", exiting..."<<endl; return;}

  PbPb = isPbPb(code);
  sample = getSample(code);
  jettree = getjettree(code);

  Init();

  outputfilenamedj = outputfolder+"/"+code+"_djt.root";
  outputfilenameinc = outputfolder+"/"+code+"_inc.root";
  outputfilenameevt = outputfolder+"/"+code+"_evt.root";



  //if (!PbPb) { do_buildtuplemc(code); return; }

  //PbPb
  //make sure folder exists!!!
  
  TString datafile = outputfolder+"dtPbj40akVs4PF_djt.root";
  if (PbPb && !file_exist(datafile)) {
    cout<<"Data file "<<outputfolder+datafile<<" doesn\'t exist. How do I calculate the centrality?"<<endl;
    return;
  }
  //make ntuples without centrality weight
  do_buildtuplemc(code);
  
  vector<float> cweight(201);
  for (int i=0;i<cweight.size();i++) cweight[i]=1;

  //build the centrality file
  if (PbPb) cweight = centrWeighting(datafile, outputfilenamedj);
  //for (int i=0;i<cweight.size();i++) cout<<i<<" - "<<cweight[i]<<endl;;
  

  update(outputfilenamedj,cweight);
  update(outputfilenameinc,cweight);
  update(outputfilenameevt,cweight);

  //repeat with good centrality file
  //do_buildtuplemc(code,centralityfile);
  
}
