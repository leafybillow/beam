#include <TObject.h>
#include <RooInt.h>
#include "TH1.h"
#include <vector>

using namespace std;
struct AHit{
  double fPosition; // Hit position in this coordinate, unit:mm
  double fCharge; // amount of charge integrated over a (isolated) cluster, unit: adc
  double fRes;  // spatial resolution of this hits, unit um
  int fMpl; // Multiplicity, unit: # of strips, a integer
};

class BeamGEMProjection: public TObject{
 private:
  vector< AHit > vHits;
  int nStrips; // should be either 256 or 512 for slac beam test
  int nHits;
  bool isSplit; // is any splitting peak ?
  TString strName; // either x or y

 public:
  BeamGEMProjection(TH1D* h, double *rms, TString name);
  ~BeamGEMProjection();

  double* fRMS; // commond mode corrected RMS array.
  TH1D* h_proj;
  
  inline vector< AHit> GetHits() const {return vHits;};
  /* inline vector<double> GetPosition() const{ return fPosition;}; */
  /* inline vector<double> GetCharge() const{ return fCharge;}; */
  /* inline vector<double> GetResolution() const {return fRes;} */
  /* inline vector<double> GetMultiplicity() const{ return fMpl;}; */
  inline bool GetSplitFlag() const {return isSplit;};
  inline int GetNHits() const {return nHits;};

  vector< pair<int,int> > SearchClusters();
  
  double ProcessCentroid(pair<int,int>);
  double ProcessCharge(pair<int,int>); 
  double ProcessResolution(pair<int,int>);
  int ProcessMultiplicity( pair<int,int>);
  

  void PlotResults();
  int CheckNStrips();

  void FitCluster(int nPeaks);
  double CheckSplit();

  void Init();
  int Process();
  int PostProcess();

  ClassDef(BeamGEMProjection,0);
};
