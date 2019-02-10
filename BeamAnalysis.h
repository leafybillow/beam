#include "TString.h"
#include "TFile.h"
#include "TH1.h"

class BeamConfig;
class BeamAnalysis{

 private:
  Int_t anaType;
  TFile *rf_raw;
  TFile *rf_output;
  
  // Tuning Parameters
  double fZSThreshold;
  
  BeamConfig *fConfig;

  Bool_t kPlot;

  int CalculatePed();
  int CalculateRMS();
  int Analysis();

  void GaussianFit(TH1D *h_fit, Double_t &mean, Double_t &sigma,
		   int iproj, int strip);

public:
  BeamAnalysis(BeamConfig *config);
  virtual ~BeamAnalysis();

  int Process();
  void PrintSummary();
  
  ClassDef(BeamAnalysis,0);
};
