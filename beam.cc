#include "BeamConfig.h"
#include "BeamAnalysis.h"

#include <iostream>

using namespace std;

void PrintUsage();

int main(int argc, char **argv){
  int run_num;
  char *configName;
  char *runType;
  char *output_name;
  char *input_name;

  bool kRunNumDefine = 0;
  bool kConfigNameDefine = 0;
  bool kRunTypeDefine = 0;
  bool kInputDefine = 0;
  bool kOutputDefine = 0;
  bool kPlot = 0;
  bool kStatus  =1;

  if(argc ==1){
    PrintUsage();
    return 0;
  }

  for(int iarg=1;iarg<argc;iarg++){
    if(strcmp(argv[iarg],"-h")==0){
      PrintUsage();
      return 0;
    }
    if(strcmp(argv[iarg],"-r")==0){
      run_num = atoi(argv[++iarg]);
      kRunNumDefine = 1;
    }
    else if(strcmp(argv[iarg],"-c")==0){
      configName = argv[++iarg];
      kConfigNameDefine =1;
    }
    else if(strcmp(argv[iarg],"-t")==0){
      runType = argv[++iarg];
      kRunTypeDefine =1;
    }
    else if(strcmp(argv[iarg],"-f")==0){
      input_name = argv[++iarg];
      kInputDefine=1;
    }
    else if(strcmp(argv[iarg],"-o")==0){
      output_name = argv[++iarg];
      kOutputDefine=1;
    }
    else if(strcmp(argv[iarg],"-P")==0){
      kPlot =1;
    }
    else{
      cerr<<__FILE__<<": "
	  <<__FUNCTION__<<": "
	  <<"unknown flag: " << argv[iarg] <<endl;
      cerr<< "See beam -h " << endl;
      kStatus =0;
      break;
    }
  }
  
  BeamConfig* fConfig = new BeamConfig();
  
  if(kConfigNameDefine)
    fConfig->SetConfigFile(Form(configName));
  else
    fConfig->SetConfigFile("beam.conf"); 

  if(!kRunNumDefine && !kInputDefine){
    cerr<< "Fatal Error: No input data specified. See beam -h."<<endl;
    kStatus=0;
  }
  else if(kRunNumDefine && kInputDefine){
    cerr<< "Fatal Error: Too many inputs defined. See beam -h."<<endl;
    kStatus=0;
  }
  else if(kRunNumDefine){
    fConfig->SetRunNumber(run_num);
  }
  else if (kInputDefine){
    fConfig->SetInputName(input_name);
  }

  if(kOutputDefine)
    fConfig->SetOutputName(output_name);

  if(kPlot)
    fConfig->SetPlotMode(kPlot);
  
  int anaType =0; // ana mode by default
  if(kRunTypeDefine){
    if(strcmp(runType,"ana")==0)
      anaType = 0;
    else if(strcmp(runType,"ped")==0)
      anaType = 1;
    else if(strcmp(runType,"rms")==0)
      anaType = 2;
    else{
      cerr<<"Error: unknown analysis type.  " <<runType << endl;
      cerr<<"See beam -h"<<endl;
      kStatus=0;
    }
  }
  
  if (kStatus)
    fConfig->SetAnalysisType(anaType);
  
  if(kStatus){
    fConfig->Config();
    BeamAnalysis *bAnalysis = new BeamAnalysis(fConfig);
    bAnalysis->Process();
    return 0;
  }
  else{
    cerr<< "Failed to configure analysis. Aborted." <<endl;
    return 1;
  }
}


void PrintUsage(){
  cout << endl;  
  cout<< "Beam Test Analysis Software " << endl;
  cout<< "author : Tao Ye <tao.ye@stonybrook.edu> " << endl;
  cout << endl;  
  cout<< "Usage: beam [-t] [-r] [-f] [-c] [-o] [-P] [-h] " << endl;
  cout<<"Options:" << endl;
  cout<< "\t" << "-h : "
      <<"Print help info " << endl;
  cout<< "\t" << "-t <analysis_type>: "
      <<" Optional, available options: ped, rms, ana. " << endl;
  cout << "\t \t : Default mode, reconstruction and  tracking analysis" << endl;
  cout << "\t \t : Generate a pedestal DB file for specific run" << endl;
  cout << "\t \t : Generate a pedestal RMS table file" << endl;
  
  cout<< "\t" << "-c <config_file>: "
      <<"optional, use beam.config by default " << endl;
  cout<< "\t" << "-r <run_num>: "
      << "mandatory, if input file is not specfied " << endl;
  cout<< "\t" << "-f <input_filename>: "
      << "mandatory, if <run_num> is not specfied " << endl;
  cout<< "\t" << "-o <output_filename>: "
      <<"optional, specify output rootfile/plot name " << endl;
  cout<< "\t" << "-P: "
      <<" output plots ONLY, no rootfile will be created. " << endl;
  cout<< endl;
}
