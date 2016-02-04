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
  
  return TString((isPbPb(code1) ? "PbPb_" : "pp_"))+getSample(code1)+"_vs_"+getSample(code2)+"_"+algo(code1);
}
