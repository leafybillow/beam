#include "TString.h"

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

class BeamConfig{

 private:
  Double_t fZSThreshold; // Zero-Suppression threshold
  TString rootfile_path;
  TString rootfile_prefix;
  TString db_template;
  TString db_path;
  Int_t n_gem;
  
  Int_t run_num;
  Bool_t kPlot;
  
  TString output_name;
  TString input_name;

  Int_t analysisType;
  TString configName;
  ifstream configFile;
  
  vector<Double_t> gem_position;
  vector<Double_t> det_position;
  
  int ParseFile();
  vector<TString>  ParseLine(TString, TString);
  void PrintSummary();

 public:
  BeamConfig();
  virtual ~BeamConfig();

  void Config();
  
  inline TString GetInputName() const { return input_name;};
  inline TString GetOutputName() const { return output_name;};
  inline TString GetDBTemplate() const { return db_template;};
  inline Int_t GetRunNumber() const { return run_num;};
  inline Double_t GetZSThreshold() const { return fZSThreshold;};
  inline Int_t GetAnalysisType() const{ return analysisType;};
  inline Bool_t GetPlotMode() const{ return kPlot;};
  
  inline void SetOutputName(TString str){output_name=str;};
  inline void SetInputName(TString str){input_name=str;};
  inline void SetAnalysisType(Int_t i){analysisType=i;};
  inline void SetRunNumber(Int_t i){run_num=i;};
  inline void SetConfigFile(TString str){configName=str;};
  inline void SetPlotMode(Bool_t flag){kPlot=flag;};
  
  ClassDef(BeamConfig,0);
};
