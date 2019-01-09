#include <TObject.h>
#include <RooInt.h>

class BeamGEMStrip: public TObject{
 private:
  double fData[6]; // container for APV25 samples
  double fAmpl_raw; // Minimum point from sample with respect to offset
  int  fT_min; //  Sample index with Minimum value 
  double fIntegral; // Integral with pedestal and common mode corrections

  // following variables are given from fittig parameters
  // Function : RC-RC shaping signal
  double fAmpl_fit; // amplitude calculated
  double fTau; //rc time constant , unit: *25ns 
  double fT_start; // starting time, unit: *25ns

 public: 
  BeamGEMStrip(double* d);
  ~BeamGEMStrip();
  inline double GetRawAmplitude() const { return fAmpl_raw;};
  inline double GetIntegral() const { return fIntegral;};
  inline int GetTMin() const { return fT_min;};

  inline double GetTau() const { return fTau;};
  inline double GetTStart() const { return fT_start;};
  inline double GetFitAmplitude() const { return fAmpl_fit;};

  // Called by Init
  void WriteSamples(double* d);
  int FindMinimum();
  double ADCSum();

  // Called by Process

  static double CRRCShaping(double*, double*);
  void FitData();
  void PlotFitResult();

  void Init();
  void Process();

  ClassDef(BeamGEMStrip,0);
};
