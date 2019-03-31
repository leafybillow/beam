#include "BeamConfig.h"
#include "BeamAnalysis.h"
#include "BeamParameters.h"

#include <iostream>

using namespace std;

void PrintUsage();

int main(int argc, char **argv){
  int run_num;
  int nEvents;
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
  bool kStatus  =0;

  if(argc==1){
    PrintUsage();
    return 0;
  }
    
  int opt;
  while( (opt=getopt(argc,argv,":c:r:t:f:o::e:hP"))!=-1){
    switch(opt){
      
    case ':':
      cout << argv[optind-1]<< " requires value. " << endl;
      kStatus = 1;
      break;
    case '?':
      cerr<<__FILE__<<": "
  	  <<__FUNCTION__<<": "
  	  <<"unknown arguments: " << optopt <<endl;
      cerr<< "See beam -h " << endl;
      kStatus = 1;
      break;
    case 'h':
      PrintUsage();
      kStatus = 1;
      break;
    case 'r':
      run_num = atoi(optarg);
      kRunNumDefine = 1;
      break;
    case 'c':
      configName=optarg;
      kConfigNameDefine=1;
      break;
    case 't':
      runType=optarg;
      kRunTypeDefine=1;
      break;
    case 'f':
      input_name=optarg;
      kInputDefine=1;
      break;
    case 'o':
      output_name=optarg;
      kOutputDefine=1;
      break;
    case 'e':
      nEvents=atoi(optarg);
      break;
    case 'P':
      kPlot =1;
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
    kStatus=1;
  }
  else if(kRunNumDefine && kInputDefine){
    cerr<< "Fatal Error: Too many inputs defined. See beam -h."<<endl;
    kStatus=1;
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

  fConfig->SetTotalEvents(nEvents);
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
      kStatus=1;
    }
  }
  
  if(kStatus==0)
    fConfig->SetAnalysisType(anaType);
  
  if(kStatus==0){
    fConfig->Config();
    BeamAnalysis *bAnalysis = new BeamAnalysis(fConfig);
    bAnalysis->Process();
  }
  else{
    cerr<< "Failed to configure analysis. Aborted." <<endl;
  }
  return kStatus;
}


void PrintUsage(){
  cout << endl;  
  cout<< "Beam Test Analysis Software " << endl;
  cout<< "author : Tao Ye <tao.ye@stonybrook.edu> " << endl;
  cout << endl;  
  cout<< "Usage: beam [-t] [-r] [-f] [-c] [-o] [-e] [-P] [-h] " << endl;
  cout<<"Options:" << endl;
  cout<< "\t" << "-h : "
      <<"Print help info " << endl;
  cout<< "\t" << "-t <analysis_type>: "
      <<" Optional, available options: ped, rms, ana. " << endl;
  cout << "\t \t ana : Default mode, reconstruction and  tracking analysis" << endl;
  cout << "\t \t ped : Generate a pedestal DB file for specific run" << endl;
  cout << "\t \t rms : Generate a pedestal RMS table file" << endl;
  
  cout<< "\t" << "-c <config_file>: "
      <<"optional, use beam.config by default " << endl;
  cout<< "\t" << "-r <run_num>: "
      << "mandatory, if input file is not specfied " << endl;
  cout<< "\t" << "-f <input_filename>: "
      << "mandatory, if <run_num> is not specfied " << endl;
  cout<< "\t" << "-o <output_filename>: "
      <<"optional, specify output rootfile/plot name " << endl;
  cout<< "\t" << "-e: "
	<<"optional, total events to be analyzed. " << endl;
  cout<< "\t" << "-P: "
      <<" output plots ONLY, no rootfile will be created. " << endl;
  cout<< endl;
}
