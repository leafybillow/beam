#include "BeamConfig.h"

#include <vector>
#include <iostream>
ClassImp(BeamConfig);
using namespace std;

BeamConfig::BeamConfig(){
  run_num= -1;
  kPlot = 0;
}

BeamConfig::~BeamConfig(){
}

void BeamConfig::Config(){
  configFile.open(configName.Data());
  ParseFile();
  TString output_suffix;
  switch(analysisType){
  case 0:
    output_suffix="_analyzed.root";
    break;
  case 1:
    output_suffix="_ped.root";
    break;
  case 2:
    output_suffix="_rms.root";
    break;
  }
    
  if(run_num!=-1){
    output_name = Form("%s%s%d%s",rootfile_path.Data(),
		       rootfile_prefix.Data(),
		       run_num,
		       output_suffix.Data());
    input_name = Form("%s%s%d.root",rootfile_path.Data(),
		      rootfile_prefix.Data(),
		      run_num);
  }
  
  PrintSummary();
}

int BeamConfig::ParseFile(){
  
  TString comment = "#";
  TString sline;
  vector<TString > vecLine;


  while(sline.ReadLine(configFile) ){
    if(sline.Contains(comment))
      continue;
    sline.ReplaceAll(" ","");
    vecLine.push_back(sline);
  }

  std::vector<TString>::iterator iter = vecLine.begin();
  vector<TString> vecStr;
  while (iter!=vecLine.end()){
    vecStr = ParseLine(*iter, "=");
    iter++;
    if(vecStr[0].Contains("db_template")){
      db_template = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("db_path")){
      db_path = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("rootfile_path")){
      rootfile_path = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("rootfile_prefix")){
      rootfile_prefix = vecStr[1];
      continue;
    }
    if(vecStr[0].Contains("zs_threshold")){
      fZSThreshold = vecStr[1].Atof();
      continue;
    }
    if(vecStr[0].Contains("n_gem")){
      n_gem = vecStr[1].Atoi();
      continue;
    }
    else{
      cerr << __FILE__ << ":"
  	   << __FUNCTION__ << ":"
  	   << "unknown config parameters "
  	   << vecStr[0]
  	   << " will be ignored. " << endl;
    }
  }

  // } // End of Line loop
  configFile.close();
  return 0;
}


vector<TString> BeamConfig::ParseLine(TString sline, TString delim){

  vector<TString> ret;
  
  while(sline.Contains(delim)){
    Int_t index = sline.Index(delim);
    TString buff = sline(0,index);
    ret.push_back(buff);
    sline.Remove(0,index+1);
  }
  
  ret.push_back(sline);
  return ret;
}
  

void BeamConfig::PrintSummary(){
  cout << "Configuration Summary:  " << endl;
  cout << "--"
       << "Config file in use: "  << configName  << endl;
  cout << "--"
       << "ROOTfile path: "  << rootfile_path  << endl;
  cout << "--"
       << "DB template: "  << db_template  << endl;
  cout << "--"
       << "DB file path: "  << db_path  << endl;
  cout << "--"
       << "Number of GEM planes: "  << n_gem  << endl;
  cout << "--"
       << "Zero Suppression threshold: "  << fZSThreshold << endl;
}
