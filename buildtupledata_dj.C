/* buildtuplemc.C - makes lightweight ntuples from HiForest files
 * change samplesfolder if necessary and subfoldernames, pthats
 * change jettree if you want to use different jet algo
 * nonoverlapping = true means events from pt,hat_bin[i-1] do not overlap 
 *                  with pt,hat_bin[i]. This reduces statistics by 10-20%, 
 *                  but allows to check if large weights corrupt errors at high pt
 */

TString samplesfolder, outputfilename, jettree;
vector<TString> subfoldernames;


void Init(TString sampleType)
{
  TString outputfolder = "/data_CMS/cms/lisniak/bjet2015/";
  samplesfolder="/data_CMS/cms/mnguyen/bJet2015/data/v2/";
  subfoldernames = {sampleType};
  outputfilename = outputfolder+"data"+sampleType+"_dj.root";
  jettree = "ak4PFJetAnalyzer/t";
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


void buildtupledata_dj(TString sampleType="pp_PFLowPt")
{
  Init(sampleType);

  vector<double> weights (1);//bins);
  weights[0] = 1;
  
  cout<<"Calculating weights"<<endl;

  for (auto w:weights) cout<<w<<"\t";
  cout<<endl;

  int totentries = 0;

  //now fill histos
  TFile *fout = new TFile(outputfilename,"recreate");
  TNtuple *nt = new TNtuple("nt","nt","hltCSV60:hltCSV80:hltPFJet60:hltPFJet80:dijet:jtpt0:jtpt1:jtphi0:jtphi1:jteta0:jteta1:discr_csvSimple0:discr_csvSimple1");
  
  for (int i=0;i<subfoldernames.size();i++) {
    auto files = list_files(TString::Format("%s/%s/",samplesfolder.Data(),subfoldernames[i].Data()));

    for (auto filename:files) {
    cout<<endl<<"Processing file "<<filename<<endl;

    TFile *f = new TFile(filename);
    TString treename = f->Get(jettree) != 0 ? jettree : "ak4PFJetAnalyzer";
    TTreeReader reader(treename,f);
    TTreeReaderValue<int> nref(reader, "nref");
    TTreeReaderArray<float> rawpt(reader, "rawpt");
    TTreeReaderArray<float> jtpt(reader, "jtpt");
    TTreeReaderArray<float> jteta(reader, "jteta");
    TTreeReaderArray<float> jtphi(reader, "jtphi");
    TTreeReaderArray<float> discr_csvSimple(reader, "discr_csvSimple");


    TTreeReader readerhlt("hltanalysis/HltTree",f);
    TTreeReaderValue<int> PFJet40(readerhlt, "HLT_AK4PFJet40_Eta5p1_v1");
    TTreeReaderValue<int> PFJet60(readerhlt, "HLT_AK4PFJet60_Eta5p1_v1");
    TTreeReaderValue<int> PFJet80(readerhlt, "HLT_AK4PFJet80_Eta5p1_v1");

    TTreeReaderValue<int> CSV60(readerhlt, "HLT_AK4PFBJetBCSV60_Eta2p1_v1");
    TTreeReaderValue<int> CSV80(readerhlt, "HLT_AK4PFBJetBCSV80_Eta2p1_v1");
    
    int nev = reader.GetEntries(true);
    totentries+=nev;
    int onep = nev/100;
    int evCounter = 0;
    TTimeStamp t0;

    while (reader.Next()) {
      readerhlt.Next();
      evCounter++;
      if (evCounter%onep==0) {
	std::cout << std::fixed;
	TTimeStamp t1; 
	cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
      }

      int ind0, ind1;
      bool foundLJ=false, foundSJ = false;
      
      for (int j=0;j<*nref;j++) {
        if (!foundLJ) //looking for the leading jet
          if (abs(jteta[j])<2) {
            ind0 = j;
            foundLJ=true;
            continue;
          }
        if (foundLJ) {
          if (abs(jteta[j])<2) { // && deltaphi>2./3.*3.142) { ----- delta-phi selection has to be made later!
            ind1 = j;
            foundSJ = true;
            break;
          }
        }
      }

      vector<float> v;
      if (foundLJ && foundSJ)
        v = {(float)*CSV60, (float)*CSV80,(float)*PFJet60,(float)*PFJet80, 1, //1 = dijet
             jtpt[ind0], jtpt[ind1],
             jtphi[ind0], jtphi[ind1],
             jteta[ind0], jteta[ind1],
             discr_csvSimple[ind0], discr_csvSimple[ind1] };
      else if (foundLJ && !foundSJ)
        v = {(float)*CSV60, (float)*CSV80,(float)*PFJet60,(float)*PFJet80, 0, //0 = monojet
             jtpt[ind0], 0,
             jtphi[ind0], 0,
             jteta[ind0], 0,
             discr_csvSimple[ind0], 0};

      if (v.size()>0)
        nt->Fill(&v[0]);

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
