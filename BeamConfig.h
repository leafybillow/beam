#include "BeamParameters.h"

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

  Int_t run_num;
  Int_t nEvents;
  Int_t evt_shift;
  Bool_t kPlot;
  
  TString output_name;
  TString input_name;

  Int_t analysisType;
  TString configName;
  ifstream configFile;

  Int_t n_gem;
  vector<Double_t> gem_position;
  
  vector<Double_t> det_position_z;
  vector<Double_t> det_position_x;
  vector<Double_t> det_position_y;

  vector<Double_t> det_width_x;
  vector<Double_t> det_width_y;
  vector<Double_t> det_thickness;
  
  vector<Int_t> qdc_channel;
  
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
  inline Int_t GetTotalEvents() const { return nEvents;};
  inline Int_t GetAnalysisType() const{ return analysisType;};
  inline Bool_t GetPlotMode() const{ return kPlot;};
  inline TString GetDBPath() const { return db_path;};
  inline TString GetRFPath() const {return rootfile_path;};
  inline vector<Int_t> GetQDCChannel() const { return qdc_channel;};
  inline vector<Double_t> GetZ_GEM() const {return gem_position;};
  inline vector<Double_t> GetZ_Det() const {return det_position_z;};
  
  inline vector<Double_t> GetWidthX_Det() const {return det_width_x;};
  inline vector<Double_t> GetWidthY_Det() const {return det_width_y;};
  inline vector<Double_t> GetPositionX_Det() const {return det_position_x;};
  inline vector<Double_t> GetPositionY_Det() const {return det_position_y;};

  inline vector<Double_t> GetThickness_Det() const {return det_thickness;};
  
  inline Int_t GetNGEMs() const {return n_gem;};
  inline Int_t GetNDets() const {return det_position_z.size() ;};
  
  inline void SetOutputName(TString str){output_name=str;};
  inline void SetInputName(TString str){input_name=str;};
  inline void SetAnalysisType(Int_t i){analysisType=i;};
  inline void SetRunNumber(Int_t i){run_num=i;};
  inline void SetConfigFile(TString str){configName=str;};
  inline void SetPlotMode(Bool_t flag){kPlot=flag;};
  inline void SetTotalEvents(int nevts){nEvents = nevts;};
  
  inline void SetEventShift(int shift){evt_shift = shift;};
  inline Int_t GetEventShift() const {return evt_shift;};
  ClassDef(BeamConfig,0);
};
