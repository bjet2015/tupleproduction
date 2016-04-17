#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TRegexp.h"
#include "TNtuple.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TTreeReader.h"
#include "TTimeStamp.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

/* buildtuplemc.C - makes lightweight ntuples from HiForest files
 * change samplesfolder if necessary and subfoldernames, pthats
 * change jettree if you want to use different jet algo
 * nonoverlapping = true means events from pt,hat_bin[i-1] do not overlap 
 *                  with pt,hat_bin[i]. This reduces statistics by 10-20%, 
 *                  but allows to check if large weights corrupt errors at high pt
 */

#include <fstream>
//#include "parsecode.h"
#include "weighting.h"
#include "mergeFCRandBJT.C"

bool PbPb = false;
bool newFlavorProcess = false;
TString samplesfolder, jettree,sample;
vector<TString> subfoldernames;
vector<int> pthats;
int Npthat = 0;
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

  if (PbPb && sample=="qcV") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"qcd30","qcd50","qcd80","qcd100","qcd120"};
    pthats = {           30,     50,     80,     100,      120};
    CS     = {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
  }
  else if (PbPb && sample=="bfV") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"bfcr30","bfcr50","bfcr80","bfcr100","bfcr120"};
    pthats = {           30,     50,     80,  100,   120};
    CS     = {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
  }
  else  if (PbPb && sample=="bjV") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"bjet30","bjet50","bjet80","bjet100","bjet120"};
    pthats =           {           30,     50,     80,  100,   120};
    CS     =           {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
    filterefficiency = {6.269e-02,7.557e-02,8.640e-02,8.908e-02,9.265e-02};
  }

  else if (PbPb && sample=="qcd") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"qcd30/puTowerExclLimitV2","qcd50/puTowerExclLimitV2","qcd80/puTowerExclLimitV2","qcd100/puTowerExclLimitV2","qcd120/puTowerExclLimitV2"};
    pthats = {           30,     50,     80,     100,      120};
    CS     = {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
  }
  else if (PbPb && sample=="bfc") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"bfcr30/puTowerExclLimitV2","bfcr50/puTowerExclLimitV2","bfcr80/puTowerExclLimitV2","bfcr100/puTowerExclLimitV2","bfcr120/puTowerExclLimitV2"};
    pthats = {           30,     50,     80,  100,   120};
    CS     = {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
  }
  else  if (PbPb && sample=="bjt") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
    subfoldernames={"bjet30/puTowerExclLimitV2","bjet50/puTowerExclLimitV2","bjet80/puTowerExclLimitV2","bjet100/puTowerExclLimitV2","bjet120/puTowerExclLimitV2"};
    pthats =           {           30,     50,     80,  100,   120};
    CS     =           {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
    filterefficiency = {6.269e-02,7.557e-02,8.640e-02,8.908e-02,9.265e-02};
  }
  // else if (PbPb && sample=="qp8") {
  //   samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia8";
  //   subfoldernames={"qcd30","qcd50","qcd80","qcd120"};
  //   pthats        ={30,50,80,120};
  //   CS = {           3.455E-02,     4.068E-03,     4.959E-04,     7.096E-05};
  // }
  // else if (PbPb && sample=="qcs") {
  //   samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pythia6";
  //   subfoldernames={"qcd30/csJets","qcd50/csJets","qcd80/csJets","qcd120/csJets"};
  //   pthats        ={30,50,80,120};
  //   CS     = {3.360e-02,3.768e-03,4.418e-04,6.159e-05};
  // } 
  else  if (PbPb && sample=="pqc") { //pyquen QCD
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pyquenColl";
    subfoldernames={"qcd100/genMatchFix_puUp"};
    pthats = {100};
    CS     = {1.511e-04};
  }
  else  if (PbPb && sample=="pfc") { //pyquen B-jet FCR
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/PbPb/pyquenColl";
    subfoldernames={"bfcr100/genMatchFix_puUp"};
    pthats = {100};
    CS     = {4.831e-07};
  }
  // else  if (!PbPb && sample=="qp8") {
  //   samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia8";
  //   subfoldernames={"qcd30","qcd50","qcd80","qcd120"};
  //   pthats = {           30,     50,     80,     120};
  //   CS = {3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05};
  // } 
  // else  if (!PbPb && sample=="bp8") {
  //   samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia8";
  //   subfoldernames={"b30","b50","b80","b120"};
  //   pthats = {           30,     50,     80,     120};  
  //   CS = {3.455E-02, 4.068E-03, 4.959E-04, 7.096E-05}; 
  //   filterefficiency = {8.202e-02,6.674e-02,9.647e-02,1.051e-01}; 
  // } 
  else  if (!PbPb && sample=="qcd") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia6";
    subfoldernames={"qcd30/constSubV1_csvV2","qcd50/constSubV1_csvV2","qcd80/constSubV1_csvV2","qcd100/constSubV1_csvV2","qcd120/constSubV1_csvV2"};
    pthats = {           30,     50,     80,  100,  120}; 
    CS     = {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
  } 
  else  if (!PbPb && sample=="bjt") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia6";
    subfoldernames={"bjet30/constSubV1_csvV2","bjet50/constSubV1_csvV2","bjet80/constSubV1_csvV2","bjet100/constSubV1_csvV2","bjet120/constSubV1_csvV2"};
    pthats =           {           30,   50,   80,   100,     120};
    CS     =           {3.360e-02,3.768e-03,4.418e-04,1.511e-04,6.159e-05};
    filterefficiency = {6.269e-02,7.557e-02,8.640e-02,8.908e-02,9.265e-02};
  }
  else  if (!PbPb && sample=="bfc") {
    samplesfolder="/data_CMS/cms/mnguyen/bJet2015/mc/pp/pythia6";
    subfoldernames={"bfcr30/constSubV1_csvV2","bfcr50/constSubV1_csvV2","bfcr80/constSubV1_csvV2","bfcr100/constSubV1_csvV2","bfcr120/constSubV1_csvV2"};
    pthats = {           30,   50,   80,   100,     120};
    CS     = {1.242e-04,1.348e-05,1.468e-06,4.831e-07,1.889e-07};
  }
  else cout<<"Unknown type of collision-mc. Use PbPb/pp-qcd/bjet"<<endl;


  Npthat = pthats.size();
}

TString getfoldername(int binnumber)
{
  return TString::Format("%s/%s/merged_HiForestAOD.root",samplesfolder.Data(),subfoldernames[binnumber].Data());
}

int getind (float pthat)
{
  if (pthat<pthats[0])
    cout<<"pthat "<<pthat<<" is less than minimum "<<pthats[0];

  for (int i=1;i<Npthat;i++) 
    if (pthat<pthats[i]) return i-1;

  return Npthat-1;
  
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

vector<TString> list_files(const char *dirname, const char *exp=".*iForest.*\\.root")
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
  int bins = Npthat;
  vector<double> weights (bins);
  weights[0] = 1;


  cout<<"Calculating weights"<<endl;
  //with distributions matching
  /*
  cout<<" pthat "<<pthats[0]<<"\t"<<flush;
  for (int i=1;i<Npthat;i++) {
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
  /*  for (int i=0;i<Npthat;i++) {
    cout<<" pthat "<<pthats[i]<<"\t"<<flush;
    double x1=0;
    if (i!=Npthat-1) {
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
  
  for (int i=0;i<Npthat;i++) {                                                                                    
    cout<<" pthat "<<pthats[i]<<"\t"<<flush;

    auto files = list_files(TString::Format("%s/%s/",samplesfolder.Data(),subfoldernames[i].Data()));

    for (auto filename:files) {
      cout<<filename;
      TFile *f = new TFile(filename);
      TH1F *hpthattemp = new TH1F("hpthattemp","hpthattemp",400,0,400);
      
      auto t = (TTree *)f->Get(jettree);

      cout<<" " <<t->GetEntries()<<endl;

      t->Project("hpthattemp","pthat");
      hpthat->Add(hpthattemp);
    }
  }

  for (int i=0;i<Npthat;i++) {
    float csi = CS[i]*filterefficiency[i];

    if (i<Npthat-1) {
      float csi_1 = CS[i+1]*filterefficiency[i+1];
      weights[i] = (csi-csi_1)/hpthat->Integral(hpthat->FindBin(pthats[i]),hpthat->FindBin(pthats[i+1])-1);
    } else {
      weights[i] = csi/hpthat->Integral(hpthat->FindBin(pthats[i]),hpthat->GetNbinsX());
      cout<<"CS : "<<csi<<" : "<<hpthat->Integral(hpthat->FindBin(pthats[i]),hpthat->GetNbinsX())<<endl;
    }
  }

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
    if (sample!="qcs")
      return (*refparton_flavorProcess)[index];
    else return -999;
  } else {
    bool gsp = (*refparton_isGSP)[index];
    if (!gsp) return 1;
    else return 6;
  }
}



void do_buildtuplemc(TString code)
{
  auto weights = calculateWeights();
  int totentries = 0;

  //put only pthat weight
  TString djvars = TString("run:lumi:event:hltCaloJet40:hltCaloJet60:hltCaloJet80:hltCSV60:hltCSV80:pthat:pthatsample:sampleEventNumber:pthatweight:bin:vz:hiHF:bProdCode:dijet:bkgLeadingJet:")+
      "subid1:refpt1:rawpt1:jtpt1:jtphi1:jteta1:discr_csvV1_1:refparton_flavorForB1:refparton_flavorProcess1:svtxm1:discr_prob1:svtxdls1:svtxpt1:svtxntrk1:nsvtx1:nselIPtrk1:"+
      "subid2:refpt2:rawpt2:jtpt2:jtphi2:jteta2:discr_csvV1_2:refparton_flavorForB2:refparton_flavorProcess2:svtxm2:discr_prob2:svtxdls2:svtxpt2:svtxntrk2:nsvtx2:nselIPtrk2:dphi21:pairCode21:"+
      "subid3:refpt3:rawpt3:jtpt3:jtphi3:jteta3:discr_csvV1_3:refparton_flavorForB3:refparton_flavorProcess3:svtxm3:discr_prob3:svtxdls3:svtxpt3:svtxntrk3:nsvtx3:nselIPtrk3:dphi31:dphi32:pairCode31:pairCode32:"+
      "SLord:subidSL:refptSL:rawptSL:jtptSL:jtphiSL:jtetaSL:discr_csvV1_SL:refparton_flavorForBSL:refparton_flavorProcessSL:svtxmSL:discr_probSL:svtxdlsSL:svtxptSL:svtxntrkSL:nsvtxSL:nselIPtrkSL:dphiSL1:pairCodeSL1:"+
      "SBord:subidSB:refptSB:rawptSB:jtptSB:jtphiSB:jtetaSB:discr_csvV1_SB:refparton_flavorForBSB:refparton_flavorProcessSB:svtxmSB:discr_probSB:svtxdlsSB:svtxptSB:svtxntrkSB:nsvtxSB:nselIPtrkSB:dphiSB1:pairCodeSB1:"+
      "Signal2ord:subidSignal2:refptSignal2:rawptSignal2:jtptSignal2:jtphiSignal2:jtetaSignal2:discr_csvV1_Signal2:refparton_flavorForBSignal2:refparton_flavorProcessSignal2:svtxmSignal2:discr_probSignal2:svtxdlsSignal2:svtxptSignal2:svtxntrkSignal2:nsvtxSignal2:nselIPtrkSignal2:dphiSignal21:pairCodeSignal21:"+
      "SignalSLord:subidSignalSL:refptSignalSL:rawptSignalSL:jtptSignalSL:jtphiSignalSL:jtetaSignalSL:discr_csvV1_SignalSL:refparton_flavorForBSignalSL:refparton_flavorProcessSignalSL:svtxmSignalSL:discr_probSignalSL:svtxdlsSignalSL:svtxptSignalSL:svtxntrkSignalSL:nsvtxSignalSL:nselIPtrkSignalSL:dphiSignalSL1:pairCodeSignalSL1";


  TString incvars = TString("run:lumi:event:hltCaloJet40:hltCaloJet60:hltCaloJet80:hltCSV60:hltCSV80:bProdCode:pthat:pthatsample:pthatweight:bin:vz:hiHF:subid:refpt:rawpt:jtpt:jtphi:jteta:discr_csvV1_:")
      +"refparton_flavorForB:isPrimary:refparton_flavorProcess:svtxm:discr_prob:svtxdls:svtxpt:svtxntrk:nsvtx:nselIPtrk";

  //now fill histos
  TFile *foutdj = new TFile(outputfilenamedj,"recreate");
  TNtuple *ntdj = new TNtuple("nt","nt",djvars);
  TFile *foutinc = new TFile(outputfilenameinc,"recreate");
  TNtuple *ntinc = new TNtuple("nt","nt",incvars);

  TFile *foutevt = new TFile(outputfilenameevt,"recreate");
  TNtuple *ntevt = new TNtuple("nt","nt","pthat:pthatbin:pthatweight:bin:vz:hiHF");

  for (int i=0;i<Npthat;i++) {

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
    TTreeReaderArray<float> * csvv1 = new TTreeReaderArray<float>(reader,sample=="bjV" || sample=="qcV" || sample=="bfV" ? "discr_csvSimple" :"discr_csvV1");//(sample=="qcs" || sample=="pqc" || sample=="pfc") ? "discr_csvV1" : "discr_csvSimple");
    //discr_csvV1 in new forests


    TTreeReaderArray<int> refparton_flavorForB(reader, "refparton_flavorForB");
    TTreeReaderArray<bool> *refparton_isGSP;

    newFlavorProcess = sample!="qcV"; //true;//sample=="qcs" || sample=="bfc" || sample == "bjt" || sample=="pqc" || sample == "pfc";//(!PbPb && (sample=="bjt" || sample=="bfc")) ||

    if (sample!="qp8" && sample!="bp8" && !newFlavorProcess)
      refparton_isGSP = new TTreeReaderArray<bool>(reader, "refparton_isGSP");
    TTreeReaderArray<int> *refparton_flavorProcess = 0;
    TTreeReaderValue<int> *bProdCode = 0;



    if (newFlavorProcess) {
      bProdCode = new TTreeReaderValue<int>(reader, "bProdCode");
      if (sample!="qcs")
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
    TTreeReaderValue<float> hiHF(readerevt, "hiHF");

    TTreeReaderValue<unsigned int> run(readerevt, "run");
    TTreeReaderValue<unsigned int> lumi(readerevt, "lumi");
    TTreeReaderValue<unsigned long long> event(readerevt, "evt");


    TString calojet40triggerv2 = !PbPb ? "HLT_AK4CaloJet40_Eta5p1ForPPRef_v1" : "HLT_HIPuAK4CaloJet40_Eta5p1_v2";
    TString calojet60triggerv2 = !PbPb ? "HLT_AK4CaloJet60_Eta5p1ForPPRef_v1" : "HLT_HIPuAK4CaloJet60_Eta5p1_v2";
    TString calojet80triggerv2 = !PbPb ? "HLT_AK4CaloJet80_Eta5p1ForPPRef_v1" : "HLT_HIPuAK4CaloJet80_Eta5p1_v2";

    // TString csv60trigger = !PbPb ? "HLT_AK4PFBJetBCSV60_Eta2p1_v1"  : "HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1";
    // TString csv80trigger = !PbPb ? "HLT_AK4PFBJetBCSV80_Eta2p1_v1"  : "HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1";
    TString csv60triggerv2 = !PbPb ? "HLT_AK4PFBJetBCSV60_Eta2p1ForPPRef_v1"  : "HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v2";
    TString csv80triggerv2 = !PbPb ? "HLT_AK4PFBJetBCSV80_Eta2p1ForPPRef_v1"  : "HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v2";

    TTreeReader readerhlt("hltanalysis/HltTree",f);
    TTreeReaderValue<int> CaloJet40(readerhlt, calojet40triggerv2);
    TTreeReaderValue<int> CaloJet60(readerhlt, calojet60triggerv2);
    TTreeReaderValue<int> CaloJet80(readerhlt, calojet80triggerv2);
    TTreeReaderValue<int> CSV60(readerhlt, csv60triggerv2);
    TTreeReaderValue<int> CSV80(readerhlt, csv80triggerv2);

    auto isSignal = [&] (int N) {return subid[N]==0;}; // && refpt[N]>20;};



    int nev = reader.GetEntries(true);
    totentries+=nev;
    int onep = nev/100;
    int evCounter = 0;
    TTimeStamp t0;

    while (reader.Next()) {
      readerevt.Next();
      readerhlt.Next();
      evCounter++;

      //if (evCounter>10*onep) break; //for fast testing

      if (evCounter%onep==0) {
  std::cout << std::fixed;
  TTimeStamp t1; 
  cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
      }

      int b = *bin;
      float centrWeight = 1;//PbPb ? centrWeights[b] : 1;

      vector<float> vevt = {*pthat, (float)i, (float)weights[getind(*pthat)], (float)*bin, *vz, *hiHF, centrWeight,
          (float)weights[getind(*pthat)]*centrWeight};
      ntevt->Fill(&vevt[0]);


      int ind1=-1, ind2=-1, indSignal2=-1, ind3=-1, indSL=-1, indSignalSL=-1, indSB=-1; //indices of leading/subleading jets in jet array
      int SLord = 0, SignalSLord = 0, Signal2ord = 0, SBord = 0;
      bool foundJ1=false, foundJ2 = false, foundSignalJ2 = false, foundJ3 = false, foundSL = false, foundSignalSL = false, foundSB = false; //found/not found yet, for convenience
      bool bkgJ1 = false; // is leading jet coming from background?

      if (abs(*vz)<15) //event-level cuts, if not passed all foundXY = false
        for (int j=0;j<*nref;j++) {
          if (abs(jteta[j])>2.0) continue; //1.5 for PU!



          //for inclusive plots, subid==0 everywhere
          if (isSignal(j)) {
            vector<float> vinc = {(float)*run, (float)*lumi, (float)*event, (float)*CaloJet40,(float)*CaloJet60,(float)*CaloJet80, (float)*CSV60, (float)*CSV80, newFlavorProcess ? (float)*(*bProdCode) : NaN, *pthat, (float)pthats[i],(float)weights[getind(*pthat)],(float)*bin, *vz,*hiHF,
              (float)subid[j], refpt[j], rawpt[j],jtpt[j], jtphi[j], jteta[j], (*csvv1)[j],
              (float)refparton_flavorForB[j], getFlavorProcess(refparton_isGSP,refparton_flavorProcess,j),
              svtxm[j],discr_prob[j],svtxdls[j],svtxpt[j],(float)svtxntrk[j],(float)nsvtx[j],(float)nselIPtrk[j]};

              ntinc->Fill(&vinc[0]);
	  }

            //if background jumped above signal (or highest signal is outside acceptance) - then it's a "bad" event
            //this condition is propagated on every clause b/c we still need to loop jets for inclusive ntuple above
            if (!foundJ1 && !isSignal(j)) bkgJ1 = true;

          	//looking for the leading jet from signal
            if (!bkgJ1 && !foundJ1 && isSignal(j)) {
              ind1 = j;
              foundJ1=true;
            }

            //ind1!=j, because 2nd jet would be filled with j of the 1st
            if (!bkgJ1 && foundJ1 && !foundJ2 && ind1!=j) {
              ind2 = j;
              foundJ2 = true;
            } 

            if (!bkgJ1 && !foundSignalJ2) Signal2ord++;
            if (!bkgJ1 && foundJ1 && !foundSignalJ2 && ind1!=j && isSignal(j)) {
              indSignal2 = j;
              foundSignalJ2 = true;
            } 

            if (!bkgJ1 && foundJ1 && foundJ2 && !foundJ3 && ind2!=j) {
              ind3 = j;
              foundJ3 = true;
            }

          	//we need ordinal number of SL jets, so counting until found
          	//indSL != SLord because some jets are not in acceptance region
            if (!bkgJ1 && !foundSL) SLord++;
          	//ind1!=j otherwise SL will be = J1
            if (!bkgJ1 && foundJ1 && !foundSL && ind1!=j && (*csvv1)[j]>0.9) {
              indSL = j;
              foundSL = true;
            }

            if (!bkgJ1 && !foundSB) SBord++;
            if (!bkgJ1 && foundJ1 && !foundSB && ind1!=j && abs(refparton_flavorForB[j])==5) {
              indSB = j;
              foundSB = true;
            }

            if (!bkgJ1 && !foundSignalSL) SignalSLord++;
            if (!bkgJ1 && foundJ1 && !foundSignalSL && ind1!=j && (*csvv1)[j]>0.9 && isSignal(j)) {
              indSignalSL = j;
              foundSignalSL = true;
            }

          }

      vector<float> vdj;

      vdj = {(float)*run, (float)*lumi, (float)*event, (float)*CaloJet40,(float)*CaloJet60,(float)*CaloJet80, (float)*CSV60, (float)*CSV80,
        *pthat,(float)pthats[i], (float)evCounter-1, (float)weights[getind(*pthat)], (float)*bin, *vz,*hiHF,
        newFlavorProcess ? (float)*(*bProdCode) : NaN,
        foundJ1 && foundJ2 ? (float)1 : (float)0, (float)bkgJ1,

        foundJ1 ? (float)subid[ind1] : NaN,
        foundJ1 ? refpt[ind1] : NaN,
        foundJ1 ? rawpt[ind1] : NaN,
        foundJ1 ? jtpt[ind1] : NaN,
        foundJ1 ? jtphi[ind1] : NaN,
        foundJ1 ? jteta[ind1] : NaN,
        foundJ1 ? (*csvv1)[ind1] : NaN,
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
        foundJ2 ? (*csvv1)[ind2] : NaN,
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
        foundJ3 ? (*csvv1)[ind3] : NaN,
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
        foundSL ? (*csvv1)[indSL] : NaN,
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
        foundSL && foundJ1 ? (float)getPairCode(refparton_flavorForB[indSL],refparton_flavorForB[ind1]) : NaN,

        foundSB ? (float)SBord : NaN,
        foundSB ? (float)subid[indSB] : NaN,
        foundSB ? refpt[indSB] : NaN,
        foundSB ? rawpt[indSB] : NaN,
        foundSB ? jtpt[indSB] : NaN,
        foundSB ? jtphi[indSB] : NaN,
        foundSB ? jteta[indSB] : NaN,
        foundSB ? (*csvv1)[indSB] : NaN,
        foundSB ? (float)refparton_flavorForB[indSB] : NaN,
        foundSB ? getFlavorProcess(refparton_isGSP,refparton_flavorProcess,indSB) : NaN,
        foundSB ? svtxm[indSB] : NaN,
        foundSB ? discr_prob[indSB] : NaN,
        foundSB ? svtxdls[indSB] : NaN, 
        foundSB ? svtxpt[indSB] : NaN,
        foundSB ? (float)svtxntrk[indSB] : NaN,
        foundSB ? (float)nsvtx[indSB] : NaN,
        foundSB ? (float)nselIPtrk[indSB] : NaN,
        foundSB && foundJ1 ? acos(cos(jtphi[indSB]-jtphi[ind1])) : NaN,
        foundSB && foundJ1 ? (float)getPairCode(refparton_flavorForB[indSB],refparton_flavorForB[ind1]) : NaN,


        foundSignalJ2 ? (float)Signal2ord : NaN,
        foundSignalJ2 ? (float)subid[indSignal2] : NaN,
        foundSignalJ2 ? refpt[indSignal2] : NaN,
        foundSignalJ2 ? rawpt[indSignal2] : NaN,
        foundSignalJ2 ? jtpt[indSignal2] : NaN,
        foundSignalJ2 ? jtphi[indSignal2] : NaN,
        foundSignalJ2 ? jteta[indSignal2] : NaN,
        foundSignalJ2 ? (*csvv1)[indSignal2] : NaN,
        foundSignalJ2 ? (float)refparton_flavorForB[indSignal2] : NaN,
        foundSignalJ2 ? getFlavorProcess(refparton_isGSP,refparton_flavorProcess,indSignal2) : NaN,
        foundSignalJ2 ? svtxm[indSignal2] : NaN,
        foundSignalJ2 ? discr_prob[indSignal2] : NaN,
        foundSignalJ2 ? svtxdls[indSignal2] : NaN, 
        foundSignalJ2 ? svtxpt[indSignal2] : NaN,
        foundSignalJ2 ? (float)svtxntrk[indSignal2] : NaN,
        foundSignalJ2 ? (float)nsvtx[indSignal2] : NaN,
        foundSignalJ2 ? (float)nselIPtrk[indSignal2] : NaN,
        foundSignalJ2 && foundJ1 ? acos(cos(jtphi[indSignal2]-jtphi[ind1])) : NaN,
        foundSignalJ2 && foundJ1 ? (float)getPairCode(refparton_flavorForB[indSignal2],refparton_flavorForB[ind1]) : NaN,

        foundSignalSL ? (float)SignalSLord : NaN,
        foundSignalSL ? (float)subid[indSignalSL] : NaN,
        foundSignalSL ? refpt[indSignalSL] : NaN,
        foundSignalSL ? rawpt[indSignalSL] : NaN,
        foundSignalSL ? jtpt[indSignalSL] : NaN,
        foundSignalSL ? jtphi[indSignalSL] : NaN,
        foundSignalSL ? jteta[indSignalSL] : NaN,
        foundSignalSL ? (*csvv1)[indSignalSL] : NaN,
        foundSignalSL ? (float)refparton_flavorForB[indSignalSL] : NaN,
        foundSignalSL ? getFlavorProcess(refparton_isGSP,refparton_flavorProcess,indSignalSL) : NaN,
        foundSignalSL ? svtxm[indSignalSL] : NaN,
        foundSignalSL ? discr_prob[indSignalSL] : NaN,
        foundSignalSL ? svtxdls[indSignalSL] : NaN, 
        foundSignalSL ? svtxpt[indSignalSL] : NaN,
        foundSignalSL ? (float)svtxntrk[indSignalSL] : NaN,
        foundSignalSL ? (float)nsvtx[indSignalSL] : NaN,
        foundSignalSL ? (float)nselIPtrk[indSignalSL] : NaN,
        foundSignalSL && foundJ1 ? acos(cos(jtphi[indSignalSL]-jtphi[ind1])) : NaN,
        foundSignalSL && foundJ1 ? (float)getPairCode(refparton_flavorForB[indSignalSL],refparton_flavorForB[ind1]) : NaN


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

void update(TString filename, TF1 * fcentrWeight, TF1 *fvertexWeight, bool bfc = false)
{
  cout<<"Updating reweighting..."<<endl;
  auto f = new TFile(filename,"update");

  auto nt = (TTree *)f->Get("nt");

  float weight,cweight,pthatweight, vz;
  float bin;
  TBranch *cw = 0, *w = 0;
  if (PbPb) {
    cw =  nt->Branch("centrWeight",&cweight);
    nt->SetBranchAddress("bin",&bin);
  }
  w =  nt->Branch("weight",&weight);
  nt->SetBranchAddress( !bfc ? "pthatweight" : "fcrweight",&pthatweight);
  nt->SetBranchAddress("vz",&vz);


  int n = nt->GetEntries();
  int onep = n/100;
  for (int i=0;i<n;i++) {
    //if (i%onep==0) cout<<i/onep<<endl;
    nt->GetEntry(i);


    cweight=PbPb ? fcentrWeight->Eval(bin)  : 1; //cweights[bin-1]
    weight=pthatweight*cweight*fvertexWeight->Eval(vz);

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
  TString jetalgo = algo(code);

  Init();

  outputfilenamedj = outputfolder+"/"+code+"_djt.root";
  outputfilenameinc = outputfolder+"/"+code+"_inc.root";
  outputfilenameevt = outputfolder+"/"+code+"_evt.root";

  //in the future pthat weight can also be moved here...


  //data sample centrality will be matched to
  //in qcd - take inclusive jets
  //in b-jet - take b-jet triggered
  TString datasample = "";
  if (!PbPb) datasample="jpf";
  else datasample = "j60";

  //!!! IMPORTANT!!! Use jet60 sample for
  //sample=="qcd" || sample=="qcV" ? "j60" : "bjt";

  //use Pu algo in data if Vs is required for MC
  //auto subalgo=sub(code);
  //if (subalgo=="Vs") 

  TString datajetalgo = jetalgo;
  //Pu is the only algorithm to handle data
  datajetalgo.ReplaceAll("Vs","Pu");
  TString datafile = outputfolder+"dt"+(PbPb?"Pb":"pp")+datasample+datajetalgo+"_djt.root";

  if (!file_exist(datafile)) {
    cout<<"Data file "<<outputfolder+datafile<<" doesn\'t exist. How do I calculate the centrality and vertex?"<<endl;
    return;
  }

  //make ntuples without centrality weight
  do_buildtuplemc(code);

  TF1 *fcentrWeight = 0;
  if (PbPb) fcentrWeight = centrWeighting(datafile, outputfilenamedj);

  auto fvertexweight = vertexWeighting(datafile, outputfilenamedj);
  
  //merge bfc sample with b-jet filtered -> bfa sample
  if (sample=="bfc" || sample=="bfV") {
    TString bjtfiledjt = outputfolder+"/mc"+(PbPb?"Pb":"pp")+ (sample=="bfc"?"bjt":"bjV") +jetalgo+"_djt.root";
    TString bfafiledjt = outputfolder+"/mc"+(PbPb?"Pb":"pp")+ (sample=="bfc"?"bfa":"baV") +jetalgo+"_djt.root";
    TString bjtfileinc = outputfolder+"/mc"+(PbPb?"Pb":"pp")+ (sample=="bfc"?"bjt":"bjV") +jetalgo+"_inc.root";
    TString bfafileinc = outputfolder+"/mc"+(PbPb?"Pb":"pp")+ (sample=="bfc"?"bfa":"baV") +jetalgo+"_inc.root";
    cout<<"Merging FCR sample "<<outputfilenamedj<<endl<<"   with BJT sample "<<bjtfiledjt<<endl;
  
    mergeFCRandBJT(outputfilenamedj,bjtfiledjt, bfafiledjt);
    update(bfafiledjt,fcentrWeight,fvertexweight, true);
  
    mergeFCRandBJT(outputfilenameinc,bjtfileinc, bfafileinc);
    update(bfafileinc,fcentrWeight,fvertexweight, true);
  }


   update(outputfilenamedj,fcentrWeight,fvertexweight);
   update(outputfilenameinc,fcentrWeight,fvertexweight);
   update(outputfilenameevt,fcentrWeight,fvertexweight);

}
