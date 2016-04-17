/* buildtupledata.C - makes lightweight ntuples from HiForest files
 * change samplesfolder if necessary and subfoldernames
 * change jettree if you want to use different jet algo
 * collision = {PbPbBJet/PbPb/pp}
*/

#include "parsecode.h"


TString jettree;
vector<TString> subfoldernames;

bool subTag = false; //IF TRUE - SUBLEDING JET MUST BE TAGGED!!!

TString outputfolder = "/data_CMS/cms/lisniak/bjet2015/";
TString samplesfolder="/data_CMS/cms/mnguyen/bJet2015/data/";


const int NaN = -999;

TTree *GetTree(TFile *f, TString treename)
{
  TTree *t = (TTree *)f->Get(treename);
  //PbPb bjet pthat120 has only unsubtracted jets!!!!! 
  //TODO: figure out
  //  if (t==0) t = (TTree *)f->Get("ak4PFJetAnalyzer/t");
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


int getEvents(TString folder, TString condition)
{
  auto files = list_files(Form("%s/%s",samplesfolder.Data(),folder.Data()));
  double x0=0;

  for (auto f:files) {
    TFile *f0 = TFile::Open(f);
    TTree *t0 = GetTree(f0,"hltanalysis/HltTree");//jettree);
    x0 += t0->GetEntries(condition);
    f0->Close();
  }

  return x0;
}

vector<double> weights;

void calculateWeights()
{
  cout<<"Calculating weights"<<endl;

  TString lowsample = subfoldernames[0];
  TString highsample = subfoldernames[1];

  int Njt80 = getEvents(highsample, "HLT_AK4PFJet80_Eta5p1_v1");
  int Njt80And60 = getEvents(highsample, "HLT_AK4PFJet80_Eta5p1_v1 && HLT_AK4PFJet60_Eta5p1_v1");
  
  weights = {(double)Njt80/(double)Njt80And60, 1};
}

vector<float> calculateWeightsBjet(TString filenamedj)
{
  cout<<"Calculating b-jet trigger weights"<<endl;

  TFile *f0 = new TFile(filenamedj);
  auto nt = (TTree *)f0->Get("nt");

  int Njt80 = nt->GetEntries("triggermatched && hltCSV80");
  int Njt80and60 = nt->GetEntries("triggermatched && hltCSV80 && hltCSV60");

  vector<float> w = {(float)Njt80/Njt80and60, 1};
  return w;
}



double getweight(TString sample, int trig60, int trig80)
{
  if (sample == subfoldernames[0] && trig60 && !trig80) return weights[0]; //pp_PFLowPt"
  if (sample == subfoldernames[1] && trig80) return 1; //pp_PFHighPt"

  return -999;
}

float matchingDistance(float jtphi1, float jteta1, float jtphi2, float jteta2)
{
  return (jteta1-jteta2)*(jteta1-jteta2)/0.6/0.6 + (jtphi1-jtphi2)*(jtphi1-jtphi2)/0.3/0.3;
}

bool matches(float dist)
{
  return dist<1.0;
} 

int triggeredLeadingJetCSV(float leadjtphi, float leadjteta, vector<Double_t> &trigpt, vector<Double_t> &trigphi, vector<Double_t> &trigeta)
{

  vector<int> btagged;
  vector<float> btaggedDist;
  for (unsigned i=0;i<trigpt.size();i++)
    for (unsigned j=i+1;j<trigpt.size();j++)
      if (trigpt[i]==trigpt[j]) btagged.push_back(i);

  for (auto ind:btagged) btaggedDist.push_back(matchingDistance(leadjtphi, leadjteta, trigphi[ind], trigeta[ind]));
  int bestmatchB = std::min_element(btaggedDist.begin(), btaggedDist.end()) - btaggedDist.begin();

  if (btaggedDist.size()>0 && matches(btaggedDist[bestmatchB])) return bestmatchB;

  return -1;
}

int triggeredLeadingJetCalo(float leadjtphi, float leadjteta, vector<Double_t> &trigpt, vector<Double_t> &trigphi, vector<Double_t> &trigeta)
{
  vector<float> caloDist;
  for (unsigned i=0;i<trigpt.size();i++) caloDist.push_back(matchingDistance(leadjtphi, leadjteta, trigphi[i], trigeta[i]));
  int bestmatchCalo = std::min_element(caloDist.begin(), caloDist.end()) - caloDist.begin();
  if (caloDist.size()>0 && matches(caloDist[bestmatchCalo])) return bestmatchCalo;

  return -1;
}


void Init(bool PbPb, TString sample)
{
  if (!PbPb && sample=="jpf") {
    subfoldernames = {"pp_PFLowPt/constSubV1_csvV2","pp_PFHighPt/constSubV1_csvV2"};
    calculateWeights();
  }
  else if (PbPb && sample=="bjt") {
    subfoldernames = {"PbPb_BJetSD/constSubV2"};
  }
  else if (PbPb && sample=="j60") {
    subfoldernames = {"PbPb_Jet60"};
    weights = {1.};
  }
  // else if (PbPb && sample=="j4_") {
  //   subfoldernames = {"PbPb_Jet40/oldProd"};
  //   weights = {1.};
  // }
  else cout<<"Don\'t know collision type: PbPb"<<PbPb<<", sample "<<sample<<endl;
}


void updatePbPbBtriggerweight(TString filename, vector<float> w)
{
  auto f = new TFile(filename,"update");

  auto nt = (TTree *)f->Get("nt");

  float csv60, csv80;
  float triggermatched;
  float weight;
  TBranch *bw;

  bw =  nt->Branch("weight",&weight);

  nt->SetBranchAddress("hltCSV60",&csv60);
  nt->SetBranchAddress("hltCSV80",&csv80);
  nt->SetBranchAddress("triggermatched",&triggermatched);
  
  int n = nt->GetEntries();
  int onep = n/100;
  for (int i=0;i<n;i++) {
    if (i%onep==0) cout<<i/onep<<endl;
    nt->GetEntry(i);


    weight = 0;
    if (triggermatched && csv80) weight = w[1];
    if (triggermatched && csv60 && !csv80) weight = w[0];

    bw->Fill();
  }

  nt->Write();
  f->Close();


}


void updateweight(TString filename)
{
  auto f = new TFile(filename,"update");

  auto nt = (TTree *)f->Get("nt");

  float prew, weight;
  TBranch *bw;

  bw =  nt->Branch("weight",&weight);
  nt->SetBranchAddress("prew",&prew);
  
  int n = nt->GetEntries();
  int onep = n/100;
  for (int i=0;i<n;i++) {
    if (i%onep==0) cout<<i/onep<<endl;
    nt->GetEntry(i);
    weight = prew;
    bw->Fill();
  }

  nt->Write();
  f->Close();


}


void buildtupledata(TString code)//(TString collision = "PbPbBJet", TString jetalgo = "akVs4PFJetAnalyzer")
{
  if (!dt(code)) { cout<<"Not data: "<<code<<", exiting..."<<endl; return;}
  
  bool PbPb = isPbPb(code);
  TString sample = getSample(code);
  jettree = getjettree(code);
  subTag = subTagging(code);

  Init(PbPb, sample);

  TString outputfilenamedj = outputfolder+"/"+code+"_djt.root";
  TString outputfilenameinc = outputfolder+"/"+code+"_inc.root";
  TString outputfilenameevt = outputfolder+"/"+code+"_evt.root";

  TString djvars = TString("run:lumi:event:prew:triggermatched:bin:vz:hiHF:hltCSV60:hltCSV80:hltCaloJet40:hltCaloJet60:hltCaloJet80:hltPFJet60:hltPFJet80:dijet:")+
      "hltCalo60jtpt:hltCalo60jtphi:hltCalo60jteta:hltCalo80jtpt:hltCalo80jtphi:hltCalo80jteta:hltCSV60jtpt:hltCSV60jtphi:hltCSV60jteta:hltCSV80jtpt:hltCSV80jtphi:hltCSV80jteta:"+
      "rawpt1:jtpt1:jtphi1:jteta1:discr_csvV1_1:svtxm1:discr_prob1:svtxdls1:svtxpt1:svtxntrk1:nsvtx1:nselIPtrk1:"+
      "rawpt2:jtpt2:jtphi2:jteta2:discr_csvV1_2:svtxm2:discr_prob2:svtxdls2:svtxpt2:svtxntrk2:nsvtx2:nselIPtrk2:dphi21:"+
      "rawpt3:jtpt3:jtphi3:jteta3:discr_csvV1_3:svtxm3:discr_prob3:svtxdls3:svtxpt3:svtxntrk3:nsvtx3:nselIPtrk3:dphi31:dphi32:"+
      "SLord:rawptSL:jtptSL:jtphiSL:jtetaSL:discr_csvV1_SL:svtxmSL:discr_probSL:svtxdlsSL:svtxptSL:svtxntrkSL:nsvtxSL:nselIPtrkSL:dphiSL1";





  for (auto w:weights) cout<<w<<"\t";
  cout<<endl;

  int totentries = 0;

  //now fill histos
  TFile *foutdj = new TFile(outputfilenamedj,"recreate");
  TNtuple *ntdj = new TNtuple("nt","ntdj",djvars);

  TFile *foutinc = new TFile(outputfilenameinc,"recreate");
  TNtuple *ntinc = new TNtuple("nt","ntinc","prew:goodevent:bin:vz:hiHF:hltCSV60:hltCSV80:hltCaloJet40:hltCaloJet60:hltCaloJet80:hltPFJet60:hltPFJet80:rawpt:jtpt:jtphi:jteta:discr_csvV1:svtxm:discr_prob:svtxdls:svtxpt:svtxntrk:nsvtx:nselIPtrk");
  TFile *foutevt = new TFile(outputfilenameevt,"recreate");
  TNtuple *ntevt = new TNtuple("nt","ntinc","prew:bin:vz:hiHF:hltCSV60:hltCSV80");
  
  for (unsigned i=0;i<subfoldernames.size();i++) {
    //get all files for unmerged forests
    auto files = list_files(TString::Format("%s/%s/",samplesfolder.Data(),subfoldernames[i].Data()));

    for (auto filename:files) {
    cout<<endl<<"Processing file "<<filename<<endl;

    TFile *f = new TFile(filename);
    TString treename = jettree;//f->Get(jettree) != 0 ? jettree : "ak3PFJetAnalyzer";
    TTreeReader reader(treename,f);
    TTreeReaderValue<int> nref(reader, "nref");
    TTreeReaderArray<float> rawpt(reader, "rawpt");
    TTreeReaderArray<float> jtpt(reader, "jtpt");
    TTreeReaderArray<float> jteta(reader, "jteta");
    TTreeReaderArray<float> jtphi(reader, "jtphi");
    TTreeReaderArray<float> discr_csvV1(reader, "discr_csvV1");

    TTreeReaderArray<float> discr_prob(reader, "discr_prob");
    TTreeReaderArray<float> svtxm(reader, "svtxm");
    TTreeReaderArray<float> svtxdls(reader, "svtxdls");
    TTreeReaderArray<float> svtxpt(reader, "svtxpt");

    TTreeReaderArray<int> svtxntrk(reader, "svtxntrk");
    TTreeReaderArray<int> nsvtx(reader, "nsvtx");
    TTreeReaderArray<int> nselIPtrk(reader, "nselIPtrk");

    TTreeReaderArray<float> *muMax=0, *muMaxTRK=0, *muMaxGBL=0;
    if (PbPb) {
      muMax = new TTreeReaderArray<float> (reader, "muMax");
      muMaxTRK = new TTreeReaderArray<float>(reader, "muMaxTRK");
      muMaxGBL = new TTreeReaderArray<float>(reader, "muMaxGBL");
    }


    //HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1 HLT_HIPuAK4CaloJet80_Eta5p1_v1

    TString calojet40trigger = !PbPb ? "HLT_AK4CaloJet40_Eta5p1_v1" : "HLT_HIPuAK4CaloJet40_Eta5p1_v1";
    TString calojet40triggerv2 = !PbPb ? "HLT_AK4CaloJet40_Eta5p1_v1" : "HLT_HIPuAK4CaloJet40_Eta5p1_v2";
    TString calojet60trigger = !PbPb ? "HLT_AK4CaloJet60_Eta5p1_v1" : "HLT_HIPuAK4CaloJet60_Eta5p1_v1";
    TString calojet80trigger = !PbPb ? "HLT_AK4CaloJet80_Eta5p1_v1" : "HLT_HIPuAK4CaloJet80_Eta5p1_v1";
    //dummy vars in PbPb case
    TString pfjet60trigger = !PbPb ? "HLT_AK4PFJet60_Eta5p1_v1" : "LumiBlock";
    TString pfjet80trigger = !PbPb ? "HLT_AK4PFJet80_Eta5p1_v1" : "LumiBlock";
    TString csv60trigger = !PbPb ? "HLT_AK4PFBJetBCSV60_Eta2p1_v1"  : "HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1";
    TString csv80trigger = !PbPb ? "HLT_AK4PFBJetBCSV80_Eta2p1_v1"  : "HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1";

    //PbPb pprimaryVertexFilter && pclusterCompatibilityFilter do nothing
    vector<TString> filterNames;
    if (PbPb) filterNames = {"pcollisionEventSelection", "HBHENoiseFilterResultRun2Loose"};
    else filterNames = {"pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "pBeamScrapingFilter"}; 

    TTreeReader readerhlt("hltanalysis/HltTree",f);
    TTreeReaderValue<int> PFJet60(readerhlt, pfjet60trigger);
    TTreeReaderValue<int> PFJet80(readerhlt, pfjet80trigger);


    TTreeReaderValue<int> CaloJet40(readerhlt, calojet40trigger);
    TTreeReaderValue<int> CaloJet40v2(readerhlt, calojet40triggerv2);
    TTreeReaderValue<int> CaloJet60(readerhlt, calojet60trigger);
    TTreeReaderValue<int> CaloJet80(readerhlt, calojet80trigger);

    TTreeReaderValue<int> CSV60(readerhlt, csv60trigger);
    TTreeReaderValue<int> CSV80(readerhlt, csv80trigger);

    TTreeReader readercsv60object("hltobject/HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v",f);
    TTreeReaderValue<vector<Double_t> > csv60pt(readercsv60object, "pt");
    TTreeReaderValue<vector<Double_t> > csv60eta(readercsv60object, "eta");
    TTreeReaderValue<vector<Double_t> > csv60phi(readercsv60object, "phi");
    
    TTreeReader readercsv80object("hltobject/HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v",f);
    TTreeReaderValue<vector<Double_t> > csv80pt(readercsv80object, "pt");
    TTreeReaderValue<vector<Double_t> > csv80eta(readercsv80object, "eta");
    TTreeReaderValue<vector<Double_t> > csv80phi(readercsv80object, "phi");
    
    TTreeReader readerCalo60object("hltobject/HLT_HIPuAK4CaloJet60_Eta5p1_v",f);
    TTreeReaderValue<vector<Double_t> > calo60pt(readerCalo60object, "pt");
    TTreeReaderValue<vector<Double_t> > calo60eta(readerCalo60object, "eta");
    TTreeReaderValue<vector<Double_t> > calo60phi(readerCalo60object, "phi");
    
    TTreeReader readerCalo80object("hltobject/HLT_HIPuAK4CaloJet80_Eta5p1_v",f);
    TTreeReaderValue<vector<Double_t> > calo80pt(readerCalo80object, "pt");
    TTreeReaderValue<vector<Double_t> > calo80eta(readerCalo80object, "eta");
    TTreeReaderValue<vector<Double_t> > calo80phi(readerCalo80object, "phi");


    TTreeReader readerevt("hiEvtAnalyzer/HiTree",f);
    TTreeReaderValue<float> vz(readerevt, "vz");
    TTreeReaderValue<int> bin(readerevt, "hiBin");
    TTreeReaderValue<float> hiHF(readerevt, "hiHF");
    
    TTreeReaderValue<unsigned int> run(readerevt, "run");
    TTreeReaderValue<unsigned int> lumi(readerevt, "lumi");
    TTreeReaderValue<unsigned long long> event(readerevt, "evt");

    TTreeReader readerskim("skimanalysis/HltTree",f);

    vector<TTreeReaderValue<int> *>filters;
    for (auto f:filterNames)
      filters.push_back(new TTreeReaderValue<int>(readerskim, f));
      
    cout<<"added filters"<<endl;
    
    int nev = reader.GetEntries(true); cout<<nev<<endl;
    totentries+=nev;
    int onep = nev/100;
    int evCounter = 0;
    TTimeStamp t0;
    
    //for testing - only 10% of data
    //while (evCounter<2*onep && reader.Next()) {
    //go full file
    while (reader.Next()) {
      readerhlt.Next();
      readerevt.Next();
      readerskim.Next();
      readercsv60object.Next();
      readercsv80object.Next();
      readerCalo60object.Next();
      readerCalo80object.Next();

      evCounter++;
      if (evCounter%onep==0) {
        std::cout << std::fixed;
        TTimeStamp t1; 
        cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
      }


      int bPFJet60 = !PbPb ? *PFJet60 : 1;
      int bPFJet80 = !PbPb ? *PFJet80 : 1;

      //int jet40 = *CaloJet40 || *CaloJet40v2;

      float weight = 1;

      if (!PbPb)
        weight = getweight(subfoldernames[i], bPFJet60, bPFJet80);


      if (PbPb && sample=="j60")
        weight = *CaloJet60;//only calojet 40

      ntevt->Fill(weight, *bin, *vz, *hiHF, *CSV60, *CSV80);

      if (weight==0) continue;

      //good event is vertex cut and noise cuts
      bool goodevent = abs(*vz)<15;
      for (auto f:filters) 
        goodevent&=*(*f);

      int ind1=-1, ind2=-1, ind3=-1, indSL=-1; //indices of leading/subleading jets in jet array
      int indTrigCSV60=-1, indTrigCSV80=-1, indTrigCalo60=-1, indTrigCalo80=-1;
      int SLord = 0;
      bool foundJ1=false, foundJ2 = false, foundJ3 = false, foundSL = false; //found/not found yet, for convenience

      bool triggermatched = false;

      if (goodevent)
        for (int j=0;j<*nref;j++) {
          //acceptance selection
          if (abs(jteta[j])>1.5) continue;
          //muon cuts
          if (PbPb) {
            if((*muMax)[j]/rawpt[j]>0.95) continue;
            if( ((*muMaxTRK)[j]-(*muMaxGBL)[j]) / ((*muMaxTRK)[j]+(*muMaxGBL)[j]) > 0.1) continue;
          }
  
          if (!foundJ1) { //looking for the leading jet
              ind1 = j;
              foundJ1=true;

	      if (PbPb) {
		indTrigCSV60 = triggeredLeadingJetCSV(jtphi[j], jteta[j], *csv60pt, *csv60phi, *csv60eta);
		indTrigCSV80 = triggeredLeadingJetCSV(jtphi[j], jteta[j], *csv80pt, *csv80phi, *csv80eta);
		indTrigCalo60 = triggeredLeadingJetCalo(jtphi[j], jteta[j], *calo60pt, *calo60phi, *calo60eta);
		indTrigCalo80 = triggeredLeadingJetCalo(jtphi[j], jteta[j], *calo80pt, *calo80phi, *calo80eta);
	      }
             
	      triggermatched = !PbPb || indTrigCSV60!=-1 || indTrigCSV80!=-1;
	  } else
            if (foundJ1 && !foundJ2) {
              ind2 = j;
              foundJ2 = true;
            } else
            if (foundJ1 && foundJ2 && !foundJ3) {
              ind3 = j;
              foundJ3 = true;
            }

          //we need ordinal number of SL jets, so counting until found
          //indSL != SLord because some jets are not in acceptance region
            if (!foundSL) SLord++;

          //ind1!=j otherwise SL will be = J1
            if (foundJ1 && ind1!=j && !foundSL && discr_csvV1[j]>0.9) {
              indSL = j;
              foundSL = true;
            }  




          //at this point foundLJ = true always, so triggermatched is determined
          vector<float> vinc = {weight, (float)triggermatched, (float) *bin, *vz, *hiHF,(float)*CSV60, (float)*CSV80,(float)*CaloJet40, (float)*CaloJet60, (float)*CaloJet80,
            (float)bPFJet60,(float)bPFJet80, rawpt[j], jtpt[j], jtphi[j], jteta[j], discr_csvV1[j],svtxm[j],discr_prob[j],
            svtxdls[j],svtxpt[j],(float)svtxntrk[j],(float)nsvtx[j],(float)nselIPtrk[j]};
  
          ntinc->Fill(&vinc[0]);
        }

      //fill dijet ntuple
      vector<float> vdj;

      vdj = {(float)*run, (float)*lumi, (float)*event, weight, (float)triggermatched, (float)*bin, *vz,*hiHF,
        (float)*CSV60, (float)*CSV80,(float)*CaloJet40,(float)*CaloJet60, (float)*CaloJet80,(float)bPFJet60,(float)bPFJet80, 
        foundJ1 && foundJ2 ? (float)1 : (float)0,

        indTrigCalo60!=-1 ? (float)(*calo60pt)[indTrigCalo60] : NaN,
        indTrigCalo60!=-1 ? (float)(*calo60phi)[indTrigCalo60] : NaN,
        indTrigCalo60!=-1 ? (float)(*calo60eta)[indTrigCalo60] : NaN,

        indTrigCalo80!=-1 ? (float)(*calo80pt)[indTrigCalo80] : NaN,
        indTrigCalo80!=-1 ? (float)(*calo80phi)[indTrigCalo80] : NaN,
        indTrigCalo80!=-1 ? (float)(*calo80eta)[indTrigCalo80] : NaN,

        indTrigCSV60!=-1  ? (float)(*csv60pt)[indTrigCSV60] : NaN,
        indTrigCSV60!=-1  ? (float)(*csv60phi)[indTrigCSV60] : NaN,
        indTrigCSV60!=-1  ? (float)(*csv60eta)[indTrigCSV60] : NaN,

        indTrigCSV80!=-1  ? (float)(*csv80pt)[indTrigCSV80] : NaN,
        indTrigCSV80!=-1  ? (float)(*csv80phi)[indTrigCSV80] : NaN,
        indTrigCSV80!=-1  ? (float)(*csv80eta)[indTrigCSV80] : NaN,
                                
        foundJ1 ? rawpt[ind1] : NaN,
        foundJ1 ? jtpt[ind1] : NaN,
        foundJ1 ? jtphi[ind1] : NaN,
        foundJ1 ? jteta[ind1] : NaN,
        foundJ1 ? discr_csvV1[ind1] : NaN,
        foundJ1 ? svtxm[ind1] : NaN,
        foundJ1 ? discr_prob[ind1] : NaN,
        foundJ1 ? svtxdls[ind1] : NaN,
        foundJ1 ? svtxpt[ind1] : NaN,
        foundJ1 ? (float)svtxntrk[ind1] : NaN,
        foundJ1 ? (float)nsvtx[ind1] : NaN,
        foundJ1 ? (float)nselIPtrk[ind1] : NaN,

        foundJ2 ? rawpt[ind2] : NaN,
        foundJ2 ? jtpt[ind2] : NaN,
        foundJ2 ? jtphi[ind2] : NaN,
        foundJ2 ? jteta[ind2] : NaN,
        foundJ2 ? discr_csvV1[ind2] : NaN,
        foundJ2 ? svtxm[ind2] : NaN,
        foundJ2 ? discr_prob[ind2] : NaN,
        foundJ2 ? svtxdls[ind2] : NaN, 
        foundJ2 ? svtxpt[ind2] : NaN,
        foundJ2 ? (float)svtxntrk[ind2] : NaN,
        foundJ2 ? (float)nsvtx[ind2] : NaN,
        foundJ2 ? (float)nselIPtrk[ind2] : NaN,
        foundJ2 && foundJ1 ? acos(cos(jtphi[ind2]-jtphi[ind1])) : NaN,
    
        foundJ3 ? rawpt[ind3] : NaN,
        foundJ3 ? jtpt[ind3] : NaN,
        foundJ3 ? jtphi[ind3] : NaN,
        foundJ3 ? jteta[ind3] : NaN,
        foundJ3 ? discr_csvV1[ind3] : NaN,
        foundJ3 ? svtxm[ind3] : NaN,
        foundJ3 ? discr_prob[ind3] : NaN,
        foundJ3 ? svtxdls[ind3] : NaN, 
        foundJ3 ? svtxpt[ind3] : NaN,
        foundJ3 ? (float)svtxntrk[ind3] : NaN,
        foundJ3 ? (float)nsvtx[ind3] : NaN,
        foundJ3 ? (float)nselIPtrk[ind3] : NaN,
        foundJ3 && foundJ1 ? acos(cos(jtphi[ind3]-jtphi[ind1])) : NaN,
        foundJ3 && foundJ2 ? acos(cos(jtphi[ind3]-jtphi[ind2])) : NaN,

        foundSL ? (float)SLord : NaN,
        foundSL ? rawpt[indSL] : NaN,
        foundSL ? jtpt[indSL] : NaN,
        foundSL ? jtphi[indSL] : NaN,
        foundSL ? jteta[indSL] : NaN,
        foundSL ? discr_csvV1[indSL] : NaN,
        foundSL ? svtxm[indSL] : NaN,
        foundSL ? discr_prob[indSL] : NaN,
        foundSL ? svtxdls[indSL] : NaN, 
        foundSL ? svtxpt[indSL] : NaN,
        foundSL ? (float)svtxntrk[indSL] : NaN,
        foundSL ? (float)nsvtx[indSL] : NaN,
        foundSL ? (float)nselIPtrk[indSL] : NaN,
        foundSL && foundJ1 ? acos(cos(jtphi[indSL]-jtphi[ind1])) : NaN};


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
  //PutInCbins(outputfolder, code, {{0,40}, {80,200}});

  if (PbPb && sample=="bjt"){
    auto w = calculateWeightsBjet(outputfilenamedj);

    updatePbPbBtriggerweight(outputfilenamedj,w);
    updatePbPbBtriggerweight(outputfilenameinc,w);
    updatePbPbBtriggerweight(outputfilenameevt,w);
  } else {
    updateweight(outputfilenamedj);
    updateweight(outputfilenameinc);
    updateweight(outputfilenameevt);
    }

}
