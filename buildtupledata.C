/* buildtupledata.C - makes lightweight ntuples from HiForest files
 * change samplesfolder if necessary and subfoldernames
 * change jettree if you want to use different jet algo
*/

TString samplesfolder, outputfilenameinc, outputfilenamedj, jettree;
vector<TString> subfoldernames;


void Init(TString sampleType, TString jetalgo)
{
  TString outputfolder = "/data_CMS/cms/lisniak/bjet2015/";
  samplesfolder="/data_CMS/cms/mnguyen/bJet2015/data/";
  subfoldernames = {sampleType};
  outputfilenamedj = outputfolder+"data"+sampleType+jetalgo+"_dj.root";
  outputfilenameinc = outputfolder+"data"+sampleType+jetalgo+"_inc.root";
  jettree = TString::Format("%s/t",jetalgo.Data());
}

TTree *GetTree(TFile *f, TString treename)
{
  TTree *t = (TTree *)f->Get(treename);
  //PbPb bjet pthat120 has only unsubtracted jets!!!!! 
  //TODO: figure out
  //  if (t==0) t = (TTree *)f->Get("ak4PFJetAnalyzer/t");
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

  return names;
}


void buildtupledata(TString sampleType="pp_PFLowPt", TString jetalgo = "ak3PFJetAnalyzer")
{
  Init(sampleType, jetalgo);

  vector<double> weights (1);//for now
  weights[0] = 1;
  
  cout<<"Calculating weights"<<endl;

  for (auto w:weights) cout<<w<<"\t";
  cout<<endl;

  int totentries = 0;

  //now fill histos
  TFile *foutdj = new TFile(outputfilenamedj,"recreate");
  TNtuple *ntdj = new TNtuple("nt","ntdj","hltCSV60:hltCSV80:hltCaloJet60:hltCaloJet80:hltPFJet60:hltPFJet80:dijet:rawpt0:rawpt1:jtpt0:jtpt1:jtphi0:jtphi1:jteta0:jteta1:discr_csvSimple0:discr_csvSimple1");
  TFile *foutinc = new TFile(outputfilenameinc,"recreate");
  TNtuple *ntinc = new TNtuple("nt","ntinc","hltCSV60:hltCSV80:hltCaloJet60:hltCaloJet80:hltPFJet60:hltPFJet80:rawpt:jtpt:jtphi:jteta:discr_csvSimple");
  
  for (int i=0;i<subfoldernames.size();i++) {
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


    TTreeReader readerhlt("hltanalysis/HltTree",f);
    TTreeReaderValue<int> PFJet60(readerhlt, "HLT_AK4PFJet60_Eta5p1_v1");
    TTreeReaderValue<int> PFJet80(readerhlt, "HLT_AK4PFJet80_Eta5p1_v1");

    TTreeReaderValue<int> CaloJet60(readerhlt, "HLT_AK4CaloJet60_Eta5p1_v1");
    TTreeReaderValue<int> CaloJet80(readerhlt, "HLT_AK4CaloJet80_Eta5p1_v1");

    TTreeReaderValue<int> CSV60(readerhlt, "HLT_AK4PFBJetBCSV60_Eta2p1_v1");
    TTreeReaderValue<int> CSV80(readerhlt, "HLT_AK4PFBJetBCSV80_Eta2p1_v1");

    TTreeReader readerevt("hiEvtAnalyzer/HiTree",f);
    TTreeReaderValue<float> vz(readerevt, "vz");

    TTreeReader readerskim("skimanalysis/HltTree",f);
    TTreeReaderValue<int> pPAprimaryVertexFilter(readerskim, "pPAprimaryVertexFilter");
    TTreeReaderValue<int> HBHENoiseFilterResultRun2Loose(readerskim, "HBHENoiseFilterResultRun2Loose");
    TTreeReaderValue<int> pBeamScrapingFilter(readerskim, "pBeamScrapingFilter");
    
    
    int nev = reader.GetEntries(true);
    totentries+=nev;
    int onep = nev/100;
    int evCounter = 0;
    TTimeStamp t0;
    
    //for testing - only 2% of data
    //    while (evCounter<2*onep && reader.Next()) {
    //go full file
    while (reader.Next()) {
      readerhlt.Next();
      readerevt.Next();
      readerskim.Next();
      evCounter++;
      if (evCounter%onep==0) {
	std::cout << std::fixed;
	TTimeStamp t1; 
	cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
      }

      //event selection
      if (! (abs(*vz)<15 && 
	     *pPAprimaryVertexFilter && 
	     *HBHENoiseFilterResultRun2Loose && 
	     *pBeamScrapingFilter )) continue;

      int ind0, ind1; //indices of leading/subleading jets in jet array
      bool foundLJ=false, foundSJ = false; //found/not found yet, for convenience


      for (int j=0;j<*nref;j++) {
        //acceptance selection
        if (abs(jteta[j])>2) continue;

        //fill inclusive jet ntuple for every jet in the acceptance region
        vector<float> vinc;
        vinc = {(float)*CSV60, (float)*CSV80,(float)*CaloJet60, (float)*CaloJet80,(float)*PFJet60,(float)*PFJet80,
		rawpt[j], jtpt[j], jtphi[j], jteta[j], discr_csvSimple[j]};
        ntinc->Fill(&vinc[0]);

        if (!foundLJ) { //looking for the leading jet
            ind0 = j;
            foundLJ=true;
            continue;
          }
        if (foundLJ && !foundSJ) {
            ind1 = j;
            foundSJ = true;
        }


      }

      //fill dijet ntuple
      vector<float> vdj;
      if (foundLJ && foundSJ)
        vdj = {(float)*CSV60, (float)*CSV80,(float)*CaloJet60, (float)*CaloJet80,(float)*PFJet60,(float)*PFJet80, 1, //1 = dijet
             rawpt[ind0], rawpt[ind1],
             jtpt[ind0], jtpt[ind1],
             jtphi[ind0], jtphi[ind1],
             jteta[ind0], jteta[ind1],
             discr_csvSimple[ind0], discr_csvSimple[ind1] };
      else if (foundLJ && !foundSJ)
        vdj = {(float)*CSV60, (float)*CSV80,(float)*CaloJet60, (float)*CaloJet80,(float)*PFJet60,(float)*PFJet80, 0, //0 = monojet
             rawpt[ind0], 0,
             jtpt[ind0], 0,
             jtphi[ind0], 0,
             jteta[ind0], 0,
             discr_csvSimple[ind0], 0};

      if (vdj.size()>0)
        ntdj->Fill(&vdj[0]);



    }

    f->Close();
    }
  }
  
  foutdj->cd();
  ntdj->Write();
  foutdj->Close();

  foutinc->cd();
  ntinc->Write();
  foutinc->Close();

  cout<<endl;
  cout<<"Total input entries "<<totentries<<endl;

}
