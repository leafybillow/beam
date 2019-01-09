#include <TObject.h>
#include <RooInt.h>
#include "TH1.h"
#include <vector>

using namespace std;

class BeamGEMProjection: public TObject{
 private:
  vector<double> fPosition; // Hit position in this coordinate, unit:mm
  vector<double> fCharge; // amount of charge integrated over a (isolated) cluster, unit: adc
  vector<double> fRes;  // spatial resolution of this hits
  vector<double> fMpl; // Multiplicity
  
  double *fRMS; // commond mode corrected RMS array.
  int nStrips; // should be either 256 or 512 for slac beam test
  int nHits;
  bool isSplit; // is a splitting peak ?
  TH1D* h_proj;

 public:
  BeamGEMProjection(TH1D* h, double *rms);
  ~BeamGEMProjection();

  inline vector<double> GetPosition() const{ return fPosition;};
  inline vector<double> GetCharge() const{ return fCharge;};
  inline vector<double> GetResolution() const {return fRes;}
  inline vector<double> GetMultiplicity() const{ return fMpl;};
  inline bool CheckSplit() const {return isSplit;};
  inline int GetNHits() const {return nHits;};

  vector<pair<int,int> > SearchClusters();
  
  double ProcessCentroid();
  double ProcessCharge();
  double ProcessResolution();
  double ProcessMultiplicity();
  
  void FitCluster();

  void Init();
  void Process();

  ClassDef(BeamGEMProjection,0);
};
