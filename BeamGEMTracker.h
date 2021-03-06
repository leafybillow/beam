#include <TObject.h>
#include <vector>
#include "BeamTypes.h"

using namespace std;
class BeamGEMPlane;
class BeamGEMProjection;
class BeamConfig;
class TLinearFitter;



class BeamGEMTracker: public TObject{
 private:
  
  vector<double> fSlope_zx;// [iTrack]
  vector<double> fSlope_zy;
  vector<double> fTheta;
  vector<double> fPhi;
  
  vector<vector<double> > fDet_x; //Extrapolated hits positions on detector plane
  vector<vector<double> > fDet_y; // [iDet][iTrack]

  vector<vector<double> > fDet_theta; //Extrapolated angle on detector plane
  vector<vector<double> > fDet_phi; // [iDet][iTrack]

  vector<double> fDet_z;
  vector<double> fDet_pos_x;
  vector<double> fDet_pos_y;
  vector<double> fDet_width_x;
  vector<double> fDet_width_y;
  vector<double> fDet_thickness;
  vector<double> qdc_value;
  
  // To-do  vector<double> rotation;
  vector<double> fGEM_z;
  vector< vector< double> > fHit_x; // [igem][ihit]
  vector< vector< double> > fHit_y;
  vector< int > effNhits; 
  vector<BeamGEMPlane* > vPlanes;
  vector< ATrack > vTracks; // vector of traces
  
  int nPlanes;
  int nTracks;
  int nPrimaries; // number of primaries tracks after combining delta tracks.
  
  bool isGoldenTrack;
  bool isFound;
  
  int track_npt; // number of planes available for tracking

  TLinearFitter* lf;
  
  void Init();
  bool FitATrack(ATrack* aTrack);
  ATrack GenerateCandidates(int* xpattern, int* ypattern);
  void SwapHits(int iplane, int ,int);
  void ProjectHits();

  void AccumulatePrimaries(ATrack aTrack);
public:
  BeamGEMTracker();
  ~BeamGEMTracker();

  inline int GetNTracks() const {return nTracks;};
  inline int GetNPrimaries() const {return nPrimaries;};
  inline vector< vector<double> > GetDetX() const {return fDet_x;};
  inline vector< vector<double> > GetDetY() const {return fDet_y;};
  inline vector< vector<double> > GetDetTheta() const {return fDet_theta;};
  inline vector< vector<double> > GetDetPhi() const {return fDet_phi;};

  inline bool IsGoldenTrack() const {return isGoldenTrack;};

  void LoadDetectorGeometry(BeamConfig* fConfig);
  inline void  LoadQDC(vector<double> qdc){ qdc_value = qdc;};
  void Process();
  void PlotResults(TString, int);
  void AddPlane(BeamGEMPlane* );

  ClassDef(BeamGEMTracker,0);
};
