#include "TString.h"

bool mc(TString code)
{
  TString s = TString(code(0,2));
  if (s=="mc") return true;
  if (s=="dt") return false;

  cout<<"Wrong name, neither dt nor mc : "<<s<<endl;
  return false;
}

bool dt(TString code)
{
  TString s = TString(code(0,2));
  if (s=="dt") return true;
  if (s=="mc") return false;

  cout<<"Wrong name, neither dt nor mc : "<<s<<endl;
  return false;
}

bool isPbPb(TString code)
{
  TString s = TString(code(2,2));
  if (s=="Pb") return true;
  if (s=="pp") return false;
  //if here - something strange submitted   
  cout<<"Wrong collision type "<<s<<"?"<<endl;
  return false;
}

TString getSample(TString code)
{
  return TString(code(4,3));
}

TString niceSample(TString sample)
{
  if (sample=="j40") return "Jet40";
  if (sample=="j4_") return "Jet40old";

  if (sample=="qcd") return "qcdPythia6";
  if (sample=="qp8") return "qcdPythia8";

  return sample;
}

TString getPythia(TString sample)
{
  if (sample=="qcd") return "Pythia 6";
  if (sample=="qp8") return "Pythia 8";
  return sample;
}

TString algo(TString code)
{
  return TString(code(7,code.Length()-7));
}

TString getjettree(TString code)
{
  return algo(code)+"JetAnalyzer/t";
}

TString sub(TString code)
{
  TString a = algo(code);
  TString s = TString(a(2,2));
  if (s=="Pu" || s=="Vs") return s;
  
  //  cout<<"Unknown subtraction"<<endl; 
  return "no";
}

TString radius(TString code)
{
  TString a = algo(code);
  TString s = TString(a(2,1));
  if (s=="3" || s=="4" || s=="5") return s;
  
  s = TString(a(4,1));
  if (s=="3" || s=="4" || s=="5") return s;

  return "?";
}

TString jettype(TString code)
{
  TString a = algo(code);
  TString s = TString(a(a.Length()-2,2));
  if (s=="PF") return s;
  
  s = TString(a(a.Length()-4,4));
  if (s=="Calo") return s;

  return "?jt";
  
}

bool checkcompatibility(TString code1, TString code2)
{
  //to check if algo was same for data and mc
  return isPbPb(code1)==isPbPb(code2) && algo(code1)==algo(code2);
  //  return sub(code1)==sub(code2) && 
  //    radius(code1)==radius(code2) && 
  //    jettype(code1)==jettype(code2);
}

TString nicepairname(TString code1, TString code2)
{
  if (!checkcompatibility(code1,code2)) return "incompatible";
  
  return TString((isPbPb(code1) ? "PbPb_" : "pp_"))+niceSample(getSample(code1))+"_vs_"+niceSample(getSample(code2))+"_"+algo(code1);
}

TString nicecentralitylabel(TString cbin)
{
  if (cbin=="") return "";
  if (cbin=="0_40") return "0-20%";
  if (cbin=="80_200") return "40-100%";
  return cbin;
}

void PutInCbins(TString outputfolder, TString code, vector<vector<int> > cbins) 
{ 
  if (!isPbPb(code)) return;
  cout<<"Puttin' on the Cbin..."<<endl;
  for (auto b:cbins) {
    int low = b[0];
    int high = b[1];

    cout<<"["<<low<<","<<high<<")"<<endl;
    
    TString folder = TString::Format("%s/cbin%d_%d",outputfolder.Data(),low,high);
    gSystem->MakeDirectory(folder); //returns -1 if exists
    
    TFile *fin, *fout;
    TTree *ntin, *ntout;
    
    for (auto suff:{"_djt","_inc","_evt"}) {
      fin = new TFile(outputfolder+"/"+code+TString(suff)+".root");
      ntin = (TTree *)fin->Get("nt");
      fout = new TFile(folder+"/"+code+TString(suff)+".root","recreate");
      ntout = ntin->CopyTree(Form("bin>=%d && bin<%d",low,high));
      ntout->Write();
      fout->Close();
      fin->Close();
    }
  }

}
