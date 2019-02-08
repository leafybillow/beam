#include "TString.h"
#include "TFile.h"
#include "TH1.h"

class BeamConfig;
class BeamAnalysis{

 private:
  Int_t anaType;
  TFile *rf_raw;
  TFile *rf_output;
  
  BeamConfig *fConfig;

  int CalculatePed();
  int CalculateRMS();
  int Analysis(Bool_t kPlot);

  void GaussianFit(TH1D *h_fit, Double_t &mean, Double_t &sigma,
		   int iproj, int strip);
  
 public:
  BeamAnalysis(BeamConfig *config);
  virtual ~BeamAnalysis();

  int Process();
  void PrintSummary();
  
  ClassDef(BeamAnalysis,0);
};
