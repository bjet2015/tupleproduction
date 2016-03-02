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

// void calculateWeightsBjet()
// {
//   cout<<"Calculating weights"<<endl;


//  //if(muMax[i]/rawpt[i]>0.95) continue;
//        // if( (muMaxTRK[i]-muMaxGBL[i]) / (muMaxTRK[i]+muMaxGBL[i]) > 0.1) continue;




//   TString lowsample = subfoldernames[0];
//   //  TString highsample = subfoldernames[1];

//   int Njt80=0, Njt80And60=0;

//   for (auto s:subfoldernames) {
//     Njt80 += getEvents(s, "HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1");
//     Njt80And60 += getEvents(s, "HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v1 && HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v1");
//   }

//   weights = {(double)Njt80/(double)Njt80And60, 1};
// }

vector<float> calculateWeightsBjet(TString filenamedj)
{
  cout<<"Calculating b-jet trigger weights"<<endl;

  TFile *f0 = new TFile(filenamedj);
  auto nt = (TTree *)f0->Get("nt");

  int Njt80 = nt->GetEntries("goodevent && hltCSV80");
  int Njt80and60 = nt->GetEntries("goodevent && hltCSV80 && hltCSV60");

  vector<float> w = {(float)Njt80/Njt80and60, 1};
  return w;
}



double getweight(TString sample, int trig60, int trig80)
{
  //  if (!trig60 && !trig80) return 0;
  if (sample == "pp_PFLowPt" && trig60) return weights[0];
  if (sample == "pp_PFHighPt" && trig80) return 1;

  return 0;
}

double getweightbjet(int csv60, int csv80)
{
  if (csv80) return weights[1];
  if (csv60 && !csv80) return weights[0];

  return 0;
}

bool matches(float jtphi1, float jteta1, float jtphi2, float jteta2)
{
  return (jteta1-jteta2)*(jteta1-jteta2)/0.6/0.6 + (jtphi1-jtphi2)*(jtphi1-jtphi2)/0.3/0.3<1.0;
} 

bool goodBtaggedEvent(float leadjtphi, float leadjteta, vector<Double_t> &trigpt, vector<Double_t> &trigphi, vector<Double_t> &trigeta)
{
  //b-tagged jets are doubled in trigger array so b-jet cannot be the only one
  if (trigpt.size()==1) return false; 

  vector<int> btagged;


  for (int i=0;i<trigpt.size();i++)
    for (int j=i+1;j<trigpt.size();j++)
      if (trigpt[i]==trigpt[j]) btagged.push_back(i);

  for (auto ind:btagged) 
    if (matches(leadjtphi, leadjteta, trigphi[ind], trigeta[ind])) return true;

  return false;
}

void Init(bool PbPb, TString sample)
{
  if (!PbPb && sample=="jpf") {
    subfoldernames = {"pp_PFLowPt","pp_PFHighPt"};
    calculateWeights();
  }
  else if (PbPb && sample=="bjt") {
    subfoldernames = {"PbPb_BJetSD/newProd"};
    //calculateWeightsBjet(); weights for bjets will be updated later
  }
  else if (PbPb && sample=="j40") {
    subfoldernames = {"PbPb_Jet40/newProd"};
    weights = {1.};
  }
  else if (PbPb && sample=="j4_") {
    subfoldernames = {"PbPb_Jet40/oldProd"};
    weights = {1.};
  }
  else cout<<"Don\'t know collision type: PbPb"<<PbPb<<", sample "<<sample<<endl;
}


void updatePbPbBtriggerweight(TString filename, vector<float> w)
{
  auto f = new TFile(filename,"update");

  auto nt = (TTree *)f->Get("nt");

  float csv60, csv80;
  float goodevent;
  float weight;
  TBranch *bw;

  bw =  nt->Branch("weight",&weight);

  nt->SetBranchAddress("hltCSV60",&csv60);
  nt->SetBranchAddress("hltCSV80",&csv80);
  nt->SetBranchAddress("goodevent",&goodevent);
  
  int n = nt->GetEntries();
  int onep = n/100;
  for (int i=0;i<n;i++) {
    if (i%onep==0) cout<<i/onep<<endl;
    nt->GetEntry(i);


    weight = 0;
    if (goodevent && csv80) weight = w[1];
    if (goodevent && csv60 && !csv80) weight = w[0];

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

  for (auto w:weights) cout<<w<<"\t";
  cout<<endl;

  int totentries = 0;

  //now fill histos
  TFile *foutdj = new TFile(outputfilenamedj,"recreate");
  TNtuple *ntdj = new TNtuple("nt","ntdj","prew:goodevent:bin:hltCSV60:hltCSV80:hltCaloJet40:hltCaloJet60:hltCaloJet80:hltPFJet60:hltPFJet80:dijet:rawpt0:rawpt1:jtpt0:jtpt1:jtphi0:jtphi1:jteta0:jteta1:discr_csvSimple0:discr_csvSimple1:svtxm0:svtxm1:discr_prob0:discr_prob1:svtxdls0:svtxdls1:svtxpt0:svtxpt1:svtxntrk0:svtxntrk1:nsvtx0:nsvtx1:nselIPtrk0:nselIPtrk1");
  TFile *foutinc = new TFile(outputfilenameinc,"recreate");
  TNtuple *ntinc = new TNtuple("nt","ntinc","prew:goodevent:bin:hltCSV60:hltCSV80:hltCaloJet40:hltCaloJet60:hltCaloJet80:hltPFJet60:hltPFJet80:rawpt:jtpt:jtphi:jteta:discr_csvSimple:svtxm:discr_prob:svtxdls:svtxpt:svtxntrk:nsvtx:nselIPtrk");
  TFile *foutevt = new TFile(outputfilenameevt,"recreate");
  TNtuple *ntevt = new TNtuple("nt","ntinc","prew:bin:hltCSV60:hltCSV80");
  
  for (int i=0;i<subfoldernames.size();i++) {
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
    TTreeReaderArray<float> discr_csvSimple(reader, "discr_csvSimple");

    TTreeReaderArray<float> discr_prob(reader, "discr_prob");
    TTreeReaderArray<float> svtxm(reader, "svtxm");
    TTreeReaderArray<float> svtxdls(reader, "svtxdls");
    TTreeReaderArray<float> svtxpt(reader, "svtxpt");

    TTreeReaderArray<int> svtxntrk(reader, "svtxntrk");
    TTreeReaderArray<int> nsvtx(reader, "nsvtx");
    TTreeReaderArray<int> nselIPtrk(reader, "nselIPtrk");

    TTreeReaderArray<float> muMax(reader, "muMax");
    TTreeReaderArray<float> muMaxTRK(reader, "muMaxTRK");
    TTreeReaderArray<float> muMaxGBL(reader, "muMaxGBL");



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
    
    // TTreeReader hltObjectTreeCSV60("hltobject/HLT_HIPuAK4CaloBJetCSV60_Eta2p1_v", f);
    // TTreeReader hltObjectTreeCSV80("hltobject/HLT_HIPuAK4CaloBJetCSV80_Eta2p1_v", f);

    // TTreeReaderValue< vector<Double_t> > pt60(*hltObjectTree, "pt");
    // TTreeReaderValue< vector<Double_t> > eta60(*hltObjectTree, "eta");
    // TTreeReaderValue< vector<Double_t> > phi80(*hltObjectTree, "phi");

    
    TTreeReader readerevt("hiEvtAnalyzer/HiTree",f);
    TTreeReaderValue<float> vz(readerevt, "vz");
    TTreeReaderValue<int> bin(readerevt, "hiBin");
    

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
    
    //for testing - only 2% of data
    //    while (evCounter<10*onep && reader.Next()) {
    //go full file
    while (reader.Next()) {
      readerhlt.Next();
      readerevt.Next();
      readerskim.Next();
      readercsv60object.Next();
      readercsv80object.Next();

      evCounter++;
      if (evCounter%onep==0) {
	std::cout << std::fixed;
	TTimeStamp t1; 
	cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
      }


      int bPFJet60 = !PbPb ? *PFJet60 : 1;
      int bPFJet80 = !PbPb ? *PFJet80 : 1;

      int jet40 = *CaloJet40 || *CaloJet40v2;

      float weight = 1;

      if (!PbPb)
	weight = getweight(subfoldernames[i], bPFJet60, bPFJet80);
      //if (PbPb && sample=="bjt")
//	weight = getweightbjet(*CSV60, *CSV80);
      if (PbPb && (sample=="j40" || sample=="j4_"))
	weight = jet40;//only calojet 40

      ntevt->Fill(weight, *bin, *CSV60, *CSV80);

      if (weight==0) continue;

      bool goodevent = abs(*vz)<15;
      for (auto f:filters) {
	goodevent&=*(*f);
      }


      if (!goodevent) continue; //TODO: check!

      int ind0, ind1; //indices of leading/subleading jets in jet array
      bool foundLJ=false, foundSJ = false; //found/not found yet, for convenience
      bool goodBtagevent = false;

      for (int j=0;j<*nref;j++) {
        //acceptance selection
        if (abs(jteta[j])>2) continue;
        //TODO: add muon cuts
        //if(muMax[i]/rawpt[i]>0.95) continue;
        // if( (muMaxTRK[i]-muMaxGBL[i]) / (muMaxTRK[i]+muMaxGBL[i]) > 0.1) continue;



        //fill inclusive jet ntuple for every jet in the acceptance region


        if (foundLJ && !foundSJ && (!subTag || (subTag && discr_csvSimple[j]>0.9))) {
            ind1 = j;
            foundSJ = true;
        }

        if (!foundLJ) { //looking for the leading jet
            ind0 = j;
            foundLJ=true;
            goodBtagevent = goodBtaggedEvent(jtphi[j], jteta[j], *csv60pt, *csv60phi, *csv60eta) ||
                            goodBtaggedEvent(jtphi[j], jteta[j], *csv80pt, *csv80phi, *csv80eta);
          }

        //at this point foundLJ = true always, so goodBtagevent is determined
        vector<float> vinc = {weight, (float)goodBtagevent, (float) *bin, (float)*CSV60, (float)*CSV80,(float)*CaloJet40, (float)*CaloJet60, (float)*CaloJet80,
          (float)bPFJet60,(float)bPFJet80, rawpt[j], jtpt[j], jtphi[j], jteta[j], discr_csvSimple[j],svtxm[j],discr_prob[j],
          svtxdls[j],svtxpt[j],(float)svtxntrk[j],(float)nsvtx[j],(float)nselIPtrk[j]};

        ntinc->Fill(&vinc[0]);


      }

      //fill dijet ntuple
      vector<float> vdj;
      if (foundLJ && foundSJ)
        vdj = {weight, (float)goodBtagevent, (float)*bin, (float)*CSV60, (float)*CSV80,(float)*CaloJet40,(float)*CaloJet60, (float)*CaloJet80,(float)bPFJet60,(float)bPFJet80, 1, //1 = dijet
	       rawpt[ind0], rawpt[ind1],
	       jtpt[ind0], jtpt[ind1],
	       jtphi[ind0], jtphi[ind1],
	       jteta[ind0], jteta[ind1],
	       discr_csvSimple[ind0], discr_csvSimple[ind1],
               svtxm[ind0],svtxm[ind1],discr_prob[ind0],discr_prob[ind1],
               svtxdls[ind0],svtxdls[ind1], svtxpt[ind0], svtxpt[ind1],
               (float)svtxntrk[ind0],(float)svtxntrk[ind1],(float)nsvtx[ind0],(float)nsvtx[ind1],
               (float)nselIPtrk[ind0],(float)nselIPtrk[ind1]};
      else if (foundLJ && !foundSJ)
        vdj = {weight, (float)goodBtagevent, (float)*bin, (float)*CSV60, (float)*CSV80,(float)*CaloJet40,(float)*CaloJet60, (float)*CaloJet80,(float)bPFJet60,(float)bPFJet80, 0, //0 = monojet
             rawpt[ind0], NaN,
             jtpt[ind0], NaN,
             jtphi[ind0], NaN,
             jteta[ind0], NaN,
               discr_csvSimple[ind0], NaN,
               svtxm[ind0],NaN,discr_prob[ind0],NaN,
               svtxdls[ind0],NaN, svtxpt[ind0], NaN,
               (float)svtxntrk[ind0],NaN,(float)nsvtx[ind0],NaN,
               (float)nselIPtrk[ind0],NaN};

      if (vdj.size()>0)
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
