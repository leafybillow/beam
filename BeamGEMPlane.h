#include <TObject.h>
#include <vector>
#include <TH1.h>
#include <TCanvas.h>

#include "BeamTypes.h"

class BeamGEMProjection;
using namespace std;
class BeamGEMPlane: public TObject{
 private:
  vector<double> fPos_x; //position in GEM coordinate, unit: mm
  vector<double> fPos_y;
  vector<double> fCharge_x; // charge amplitude, unit: adc
  vector<double> fCharge_y;
  vector<double> fWidth_x; // charge amplitude, unit: adc
  vector<double> fWidth_y;
  vector<double> fPeakHeight_x; // Peak Amplitude, unit: adc
  vector<double> fPeakHeight_y;

  vector<int> fSplit_x; // Split Level
  vector<int> fSplit_y;

  vector< correlator > vCorrelator;
  
  vector<int> vHitsMask_x;  // 1 : accept , 0 : reject 
  vector<int> vHitsMask_y;
  double z_position; // position in z; Not used now
  int nHits;  //Number of hits found
  TString strPlaneName;
  BeamGEMProjection* bgProjX;
  BeamGEMProjection* bgProjY;
  
  vector< AHit > xHits;
  vector< AHit > yHits;
  
  int my_id;
  //Process functions

  int Reconstruct();  // return number of hits reconstructed

  void CollectResults();

  vector<int> GenerateKeys(int, int, int, vector<int>);
  void EvalCorrelation(correlator &aCorrelator);
  correlator GenerateCorrelator(int key1, int key2);
  void UpdateCorrelator(correlator &aCorrelator);
  
  //Init Check
  int CheckProjections();

 public:
  // Called by Users
  BeamGEMPlane();
  ~BeamGEMPlane();
  
  inline vector<double> GetPositionX() const {return fPos_x;};
  inline vector<double> GetPositionY() const {return fPos_y;};
  inline vector<double> GetChargeX() const {return fCharge_x;};
  inline vector<double> GetChargeY() const {return fCharge_y;};
  inline vector<double> GetWidthX() const {return fWidth_x;};
  inline vector<double> GetWidthY() const {return fWidth_y;};
  inline vector<double> GetPeakHeightX() const {return fPeakHeight_x;};
  inline vector<double> GetPeakHeightY() const {return fPeakHeight_y;};

  inline vector<int> GetSplitX() const {return fSplit_x;};
  inline vector<int> GetSplitY() const {return fSplit_y;};

  inline int GetNHits() const {return nHits;};

  inline BeamGEMProjection* GetProjectionX() const {return bgProjX;};
  inline BeamGEMProjection* GetProjectionY() const {return bgProjY;};
  inline double GetPositionZ() const { return z_position;};
  
  inline void SetPositionZ(double z) { z_position=z;};
  inline void SetID(int id) { my_id = id;};
  
  void AddProjectionX( BeamGEMProjection* );
  void AddProjectionY( BeamGEMProjection* );

  int Process();
  void PlotResults(TString, int);
  void PrintSummary();
  
  ClassDef(BeamGEMPlane,0);
};
