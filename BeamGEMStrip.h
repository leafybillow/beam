#include <TObject.h>

class TH1D;
class BeamGEMStrip: public TObject{
 private:
  double fData[6]; // container for APV25 samples
  double fAmpl_raw; // Minimum point from sample with respect to offset
  int  fT_max; //  Sample index with peak value 
  double fADCsum; // Integral with pedestal and common mode corrections

  bool kZeroSuppression; // if 1, it is suppressed
  // following variables are given from fittig parameters
  // Function : CR-RC shaping 
  double fAmpl_fit; // amplitude calculated
  double fTau; //rc time constant , unit: *25ns 
  double fT_start; // starting time, unit: *25ns
  double Chi2; // FIXME: to-do: chi-square 

  // Identification
  int id_strip;  // strip id: the physical position in readout plane
  TH1D *h_fit;

  // Called by Init
  void WriteSamples(double* d);
  int FindMaximum();
  double SumADC();

  // Called by Process

  static double CRRCShaping(double*, double*);
  void FitData();
  void Init();

 public: 
  BeamGEMStrip(double* d, int id);
  ~BeamGEMStrip();

  inline int GetStripID() const {return id_strip;};
  inline double GetRawAmplitude() const { return fAmpl_raw;};
  inline double GetADCsum() const { return fADCsum;};
  inline int GetTMax() const { return fT_max;};
  
  inline double GetTau() const { return fTau;};
  inline double GetTStart() const { return fT_start;};
  inline double GetFitAmplitude() const { return fAmpl_fit;};

  inline double GetZSStatus() const {return kZeroSuppression;};
  
  double GetAmplitude();

  // Called by Users
  void Process();
  void PlotFitResult(int ievt,int igem, int iapv, int ich);
  void SetZeroSuppression(bool zsflag);
  
  ClassDef(BeamGEMStrip,0);
};
