#include <TObject.h>
#include <RooInt.h>
#include <vector>
using namespace std;
class BeamGEMPlane;
class BeamGEMTracker: public TObject{
 private:
  vector<double> fSlope_zx;
  vector<double> fSlope_zy;
  vector<double> fTheta;
  vector<double> fPhi;
  vector<double> fDet_x;  //Extrapolated hits positions on detector plane
  vector<double> fDet_y;
  int nTracks;
 public:
  BeamGEMTracker(BeamGEMPlane*, BeamGEMPlane*);
  ~BeamGEMTracker();

  /* void SetPositionX(vector<double> pos); */
  /* void SetPositionY(vector<double> pos); */
  /* void SetNHits(int n); */
  inline int GetNTracks() const {return nTracks;};
  
  
  void Init();
  void Process();
  ClassDef(BeamGEMTracker,0);
};
