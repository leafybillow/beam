#include <TObject.h>
#include <vector>
#include <TH1.h>
#include <TCanvas.h>

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

  vector<int> fSplit_x; // Split Level
  vector<int> fSplit_y;

  vector<double> fCorelation; // (x-y)/(x+y)

  vector<int> vHitsMask_x;  // 1 : accept , 0 : reject 
  vector<int> vHitsMask_y;
  double z_position; // position in z; Not used now
  int nHits;  //Number of hits found
  TString strPlaneName;
  BeamGEMProjection* bgProjX;
  BeamGEMProjection* bgProjY;

  //Process functions

  int Reconstruct(); // return number of hits

  void CollectResults();
  double EvalCorrelation(double charge_x, double charge_y);

  //Init Check
  int CheckProjections();
  // Post-check
  int CheckHits();  

 public:
  // Called by Users
  BeamGEMPlane(TString name);
  ~BeamGEMPlane();
  
  inline vector<double> GetPositionX() const {return fPos_x;};
  inline vector<double> GetPositionY() const {return fPos_y;};
  inline vector<double> GetChargeX() const {return fCharge_x;};
  inline vector<double> GetChargeY() const {return fCharge_y;};
  inline vector<double> GetWidthX() const {return fWidth_x;};
  inline vector<double> GetWidthY() const {return fWidth_y;};

  inline vector<int> GetSplitX() const {return fSplit_x;};
  inline vector<int> GetSplitY() const {return fSplit_y;};

  inline vector<double> GetCorrelation() const {return fCorelation;};
  inline int GetNHits() const {return nHits;};

  inline BeamGEMProjection* GetProjectionX() const {return bgProjX;};
  inline BeamGEMProjection* GetProjectionY() const {return bgProjY;};

  void AddProjectionX( BeamGEMProjection* );
  void AddProjectionY( BeamGEMProjection* );

  int Process();
  void PlotResults(TString, int);

  ClassDef(BeamGEMPlane,0);
};
