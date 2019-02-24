#include "BeamConfig.h"

#include <vector>
#include <iostream>
ClassImp(BeamConfig);
using namespace std;

// if not define, these parameters will be default
double zs_threshold = 3.0;
double width_cut=1.0;
double width_threshold = 3.0; 
double split_frac=0.1;
int edge_cut = 5;
int stability = 50;

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
      zs_threshold = vecStr[1].Atof();
      fZSThreshold = zs_threshold;
      continue;
    }
    if(vecStr[0].Contains("width_cut")){
      width_cut = vecStr[1].Atof();
      continue;
    }
    if(vecStr[0].Contains("edge_cut")){
      edge_cut = vecStr[1].Atoi();
      continue;
    }
    if(vecStr[0].Contains("stability")){
      stability = vecStr[1].Atoi();
      continue;
    }
    if(vecStr[0].Contains("width_threshold")){
      width_threshold = vecStr[1].Atof();
      continue;
    }
    if(vecStr[0].Contains("split_frac")){
      split_frac = vecStr[1].Atof();
      continue;
    }

    if(vecStr[0].Contains("n_gem")){
      n_gem = vecStr[1].Atoi();
      continue;
    }
    if(vecStr[0].Contains("gem_position")){
      vector<TString> buff = ParseLine(vecStr[1],",");
      vector<TString>::iterator iter = buff.begin();
      while(iter!=buff.end()){
	gem_position.push_back( (*iter).Atof() );
	iter++;
      }
      continue;
    }
    if(vecStr[0].Contains("det_position")){
      vector<TString> buff = ParseLine(vecStr[1],",");
      vector<TString>::iterator iter = buff.begin();
      while(iter!=buff.end()){
	det_position.push_back( (*iter).Atof() );
	iter++;
      }
      continue;
    }
    if(vecStr[0].Contains("qdc_channel")){
      vector<TString> buff = ParseLine(vecStr[1],",");
      vector<TString>::iterator iter = buff.begin();
      while(iter!=buff.end()){
	qdc_channel.push_back( (*iter).Atoi());
	iter++;
      }
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
       << "GEM Plane Z positions (mm): ";
  vector<Double_t>::iterator it_gem = gem_position.begin();
  while(it_gem!=gem_position.end() ){
    cout << *it_gem << " ";
    it_gem++;
  }
  cout << endl;
  
  cout << "--"
       << "Detector Plane Z positions (mm): ";
  vector<Double_t>::iterator it_det = det_position.begin();
  while(it_det!=det_position.end() ){
    cout << *it_det << " ";
    it_det++;
  }
  cout << endl;
  cout << "-- Analysis Parameters" << endl; 
  cout << "--"
       << "Zero Suppression threshold: "  << zs_threshold << endl;
  cout << "--"
       << "Reject cluster not greater than: "  << width_cut << endl;
  cout << "--"
       << "Cut off Number of strips in the edge : "  << edge_cut << endl;
  cout << "--"
       << "Oversize Cluster threshold: "  << width_threshold << endl;
  cout << "--"
       << "Spliting Fraction threshold: "  << split_frac << endl;
  cout << "--"
       << "Stability (adc counts): "  << stability << endl;
}
