#include <TObject.h>
#include <RooInt.h>

class BeamGEMStrip: public TObject{
 private:
  double data[6]; // container for APV25 samples
  double ampl_raw; // Minimum point from sample with respect to offset
  int  t_min;
  double integral; // Integral with pedestal and common mode corrections

  // following variables are given from fittig parameters
  // Function : RC-RC shaping signal
  double ampl_fit; // amplitude calculated
  double tau; //rc time constant , unit: *25ns 
  double t_start; // starting time, unit: *25ns

 public: 
  BeamGEMStrip(double* d);
  ~BeamGEMStrip();
  inline double GetRawAmplitude() const { return ampl_raw;};
  inline double GetIntegral() const { return integral;};
  inline int GetTMin() const { return t_min;};

  inline double GetTau() const { return tau;};
  inline double GetTStart() const { return t_start;};
  inline double GetFitAmplitude() const { return ampl_fit;};
  
  inline void SetRawAmplitude(double ampl){ ampl_raw=ampl;};
  inline void SetIntegral(double intg){ integral=intg;};
  inline void SetTMin(int t){t_min =t;};

  inline void SetFitAmplitude(double ampl){ ampl_fit =ampl;};
  inline void SetTStart(double t){t_start=t;};
  inline void SetTau(double t){tau=t;};

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
