#include <TObject.h>
#include <vector>


using namespace std;
class BeamGEMPlane;
class BeamGEMProjection;
class TLinearFitter;

struct ATrack{
  double fSlope;
  double fIntercept;
  double fChi2;
};

class BeamGEMTracker: public TObject{
 private:
  
  vector<double> fSlope_zx;// [iTrack]
  vector<double> fSlope_zy;
  vector<double> fTheta;
  vector<double> fPhi;
  
  vector<double> fDet_x; //Extrapolated hits positions on detector plane
  vector<double> fDet_y; // [iDet]
  vector<double> fDet_z;

  vector<ATrack> vTrack_zx;
  vector<ATrack> vTrack_zy;
  
  vector<double> fGEM_z;
  vector< vector< double> > fHit_x; // [igem][ihit]
  vector< vector< double> > fHit_y;
  
  vector<BeamGEMPlane* > vPlanes;
  int nPlanes;
  int nTracks;
  bool isGoldenTrack;

  int track_flag;

  TLinearFitter* lf;
  
  void Init();
  bool FitSingleTrack(int iHit);
  
 public:
  BeamGEMTracker();
  ~BeamGEMTracker();

  inline int GetNTracks() const {return nTracks;};
  
  void Process();
  void PlotResults(TString, int);
  void AddPlane(BeamGEMPlane* );

  ClassDef(BeamGEMTracker,0);
};
