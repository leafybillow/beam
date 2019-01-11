#include <TObject.h>
#include <RooInt.h>
#include <vector>
#include <TH1.h>
#include <TCanvas.h>


using namespace std;
class BeamGEMProjection;
class BeamGEMPlane: public TObject{
 private:
  vector<double> fPos_x; //position in GEM coordinate, unit: mm
  vector<double> fPos_y;
  vector<double> fCharge_x; // charge amplitude, unit: adc
  vector<double> fCharge_y;
  vector<double> fCorelation; // (x-y)/(x+y)

  double z_position; // position in z; Not used now
  int nHits;  //Number of hits found
  TString strPlaneName;
  BeamGEMProjection* bgProjX;
  BeamGEMProjection* bgProjY;
  //Process functions
  void MatchHits();
  double EvalCorrelation(double charge_x, double charge_y);
  
  //Init Check
  int CheckProjections();
  
  // Post-check
  int CheckHits();  

 public:
  BeamGEMPlane(TString);
  ~BeamGEMPlane();
  
  inline vector<double> GetPositionX() const {return fPos_x;};
  inline vector<double> GetPositionY() const {return fPos_y;};
  inline vector<double> GetChargeX() const {return fCharge_x;};
  inline vector<double> GetChargeY() const {return fCharge_y;};
  inline vector<double> GetCorrelation() const {return fCorelation;};
  inline int GetNHits() const {return nHits;};

  // Called by Users
  void AddProjectionX( BeamGEMProjection* );
  void AddProjectionY( BeamGEMProjection* );

  int Process();
  void PlotResults(TString, int);

  ClassDef(BeamGEMPlane,0);
};
