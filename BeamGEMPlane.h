#include <TObject.h>
#include <RooInt.h>
#include <vector>
#include <TH1.h>



using namespace std;
class BeamGEMProjection;
class BeamGEMPlane: public TObject{
 private:
  vector<double> fPos_x; //position in GEM coordinate, unit: mm
  vector<double> fPos_y;
  vector<double> fCharge_x; // charge amplitude, unit: adc
  vector<double> fCharge_y;
  vector<double> fCorelation; // (x-y)/(x+y)

  int nHits;  //Number of hits found
  BeamGEMProjection* gProjX;
  BeamGEMProjection* gProjY;
 public:
  BeamGEMPlane(BeamGEMProjection* projX , BeamGEMProjection* projY );
  ~BeamGEMPlane();
  
  inline vector<double> GetPositionX() const {return fPos_x;};
  inline vector<double> GetPositionY() const {return fPos_y;};
  inline vector<double> GetChargeX() const {return fCharge_x;};
  inline vector<double> GetChargeY() const {return fCharge_y;};
  inline vector<double> GetCorrelation() const {return fCorelation;};
  inline int GetNHits() const {return nHits;};

  //Process functions
  void MatchHits();
  double EvalCorrelation(double charge_x, double charge_y);
  
  //Init Check
  int CheckPlaneName();
  
  // Post-check
  int CheckHits();  

  void Init();
  int Process();

  ClassDef(BeamGEMPlane,0);
};
