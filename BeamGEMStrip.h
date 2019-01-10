#include <TObject.h>
#include <RooInt.h>
class TH1D;
class BeamGEMStrip: public TObject{
 private:
  double fData[6]; // container for APV25 samples
  double fAmpl_raw; // Minimum point from sample with respect to offset
  int  fT_max; //  Sample index with peak value 
  double fADCsum; // Integral with pedestal and common mode corrections

  // following variables are given from fittig parameters
  // Function : CR-RC shaping 
  double fAmpl_fit; // amplitude calculated
  double fTau; //rc time constant , unit: *25ns 
  double fT_start; // starting time, unit: *25ns
  double Chi2; // FIXME: to-do: chi-square 

  TH1D *h_fit;
 public: 
  BeamGEMStrip(double* d);
  ~BeamGEMStrip();
  inline double GetRawAmplitude() const { return fAmpl_raw;};
  inline double GetADCsum() const { return fADCsum;};
  inline int GetTMax() const { return fT_max;};

  inline double GetTau() const { return fTau;};
  inline double GetTStart() const { return fT_start;};
  inline double GetFitAmplitude() const { return fAmpl_fit;};

  // Called by Init
  void WriteSamples(double* d);
  int FindMaximum();
  double SumADC();

  // Called by Process

  static double CRRCShaping(double*, double*);
  void FitData();

  // Called by Users
  void Init();
  void Process();
  void PlotFitResult(int ievt,int igem, int iapv, int ich);

  ClassDef(BeamGEMStrip,0);
};
