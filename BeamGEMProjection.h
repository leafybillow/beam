#include <TObject.h>
#include <RooInt.h>
#include "TH1.h"
#include <vector>

using namespace std;
struct AHit{
  double fPosition; // Hit position in this coordinate, unit:mm
  double fCharge; // amount of charge integrated over a (isolated) cluster, unit: adc
  double fRes;  // spatial resolution of this hits, unit um
  int fWidth; // A single hit width, unit: # of strips, a integer
};

class BeamGEMStrip;

class BeamGEMProjection: public TObject{
 private:
  vector <BeamGEMStrip* > vBGStrips; // Not Used for Now
  vector< AHit > vHits;
  int nStrips; // should be either 256 or 512 for GEMs in slac beam test
  int nHits;
  bool isSplit; // is any splitting peak ?
  TString strProjName; // e.g. x1, x2 , y1, y2

  // Called by Process
  vector< pair<int,int> > SearchClusters();
  void SortHits(); 
  double ProcessCentroid(pair<int,int>);
  double ProcessCharge(pair<int,int>); 
  double ProcessResolution(pair<int,int>);
  int ProcessWidth( pair<int,int>);

  int CheckNStrips();
  double CheckSplit();

  void FitCluster(int nPeaks);
  TH1D* h_proj;

 public:
  BeamGEMProjection(TString projName, Int_t nch);
  ~BeamGEMProjection();
  
  inline vector< AHit> GetHits() const {return vHits;};
  inline bool GetSplitFlag() const {return isSplit;};
  inline int GetNHits() const {return nHits;};
  inline TString GetProjName() const{return strProjName;};
  inline int GetNStrips() const {return nStrips;};
  inline TH1D* GetTH1D() const {return h_proj;};
  void Init();
  // Called by User
  int Process();
  void AddStrip(BeamGEMStrip* );
  void PlotResults(TString, int);

  int PostProcess();
  int TestCrossTalk(int iHit1, int iHit2); // 1 : suspected as a cross talk pair; if 0: it is not 
  void UpdateHits( vector<int> vHitsMask ); // called by BeamGEMPlane::Process(), during post-check on cross talk

  ClassDef(BeamGEMProjection,0);
};
