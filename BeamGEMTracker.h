#include <TObject.h>
#include <vector>


using namespace std;
class BeamGEMPlane;
class BeamGEMProjection;
class TLinearFitter;

struct ATrace{
  vector<double> x;
  vector<double> y;
  vector<double> z;
  double fSlope_zx;
  double fSlope_zy;
  double fIntercept_x;
  double fIntercept_y;
  double fChi2;
  int myID;
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
  
  vector<double> fGEM_z;
  vector< vector< double> > fHit_x; // [igem][ihit]
  vector< vector< double> > fHit_y;
  
  vector<BeamGEMPlane* > vPlanes;
  vector< ATrace > vTraces; // vector of traces
  
  int nPlanes;
  int nTracks;
  
  bool isGoldenTrack;
  bool isFound;
  
  int track_npt;

  TLinearFitter* lf;
  
  void Init();
  bool FitATrack(ATrace* aTrace);
  void GenerateCandidates();
 public:
  BeamGEMTracker();
  ~BeamGEMTracker();

  inline int GetNTracks() const {return nTracks;};
  
  void Process();
  void PlotResults(TString, int);
  void AddPlane(BeamGEMPlane* );

  ClassDef(BeamGEMTracker,0);
};
