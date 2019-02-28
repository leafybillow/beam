#include <TObject.h>
#include "TH1.h"
#include <vector>

const int lower_neighbor[128]=
  {127,-1,120,121,122,123,124,125,
   96,97,98,99,100,101,102,103,
   104,105,106,107,108,109,110,111,
   112,113,114,115,116,117,118,119,
   0,1,2,3,4,5,6,7,
   8,9,10,11,12,13,14,15,
   16,17,18,19,20,21,22,23,
   24,25,26,27,28,29,30,31,
   32,33,34,35,36,37,38,39,
   40,41,42,43,44,45,46,47,
   48,49,50,51,52,53,54,55,
   56,57,58,59,60,61,62,63,
   64,65,66,67,68,69,70,71,
   72,73,74,75,76,77,78,79,
   80,81,82,83,84,85,86,87,
   88,89,90,91,92,93,94,95}; // [istrip]

const int higher_neighbor[128]=
  {32,33,34,35,36,37,38,39,
   40,41,42,43,44,45,46,47,
   48,49,50,51,52,53,54,55,
   56,57,58,59,60,61,62,63,
   64,65,66,67,68,69,70,71,
   72,73,74,75,76,77,78,79,
   80,81,82,83,84,85,86,87,
   88,89,90,91,92,93,94,95,
   96,97,98,99,100,101,102,103,
   104,105,106,107,108,109,110,111,
   112,113,114,115,116,117,118,119,
   120,121,122,123,124,125,126,127,
   8,9,10,11,12,13,14,15,
   16,17,18,19,20,21,22,23,
   24,25,26,27,28,29,30,31,
   2,3,4,5,6,7,-1,0}; // [istrip]
// returns physical strip number of its adjacent channels

using namespace std;

struct ACluster{
  double fPosition; // Cluster position in this projection, unit:mm
  double fCharge; // amount of charge integrated over a (isolated) cluster, unit: adc
  int fWidth; // A cluster width, unit: # of strips, an integer
  pair<int, int> pRange;
  int fSplit;  // if 0, no split is detected
  vector<int> peak; // peak position;
};

struct AHit{
  double fPosition; // Hit position in this projection, unit:mm
  double fCharge; // amount of charge integrated over this hit, unit: adc
  double fRes;  // spatial resolution of this hit, unit um
  int fWidth; // A single hit width, unit: # of strips, an integer
};
  
class BeamGEMStrip;

class BeamGEMProjection: public TObject{
 private:
  vector <BeamGEMStrip* > vBGStrips; // Not Used for Now
  
  vector< AHit > vHits;
  vector < ACluster > vClusters;
  Double_t charge_sum;
  int nStrips; // should be either 256 or 512 for GEMs in slac beam test
  int nHits;
  int nClusters;
  TString strProjName; // e.g. x1, x2 , y1, y2
  TH1D* h_proj; // zero-suppressed histogram for calculation,
  TH1D* h_raw;  // histogram for plotting results
  // use histogram for quick check
  // for fast performance, a std container is prefered
  vector<double> vStat ; // vector to collect all samples
  vector<double> vBaseline; // vector to collect baseline samples

  double baseline_mean;
  double baseline_rms;
  double overall_mean;
  double overall_rms;
  
  void Init();
  double CalculateMean(vector<double> );
  double CalculateRMS(vector<double> );
  // Called by Process
  // Coarse Process for Clusters
  int CoarseProcess();
  vector< pair<int,int> > SearchClusters();
  void SortClusters();
  double ProcessCentroid(pair<int,int>);
  double ProcessCharge(pair<int,int>); 
  int ProcessWidth( pair<int,int>);
  vector<int> ProcessSplitCheck(pair<int,int>);
  
  void RejectCrossTalk(); 
  int TestCrossTalk(ACluster i, ACluster j); // 1 : suspected as a cross talk pair; if 0: it is not
  // int TestCrossTalk_v1(int iHit1, int iHit2);

  // Fine Process for Hits
  int FineProcess();
  void SeparateHits(pair<int,int>); // FIXME: to-do
  double ProcessResolution(pair<int,int>);
  void SortHits();

  int CheckNStrips();
  void FitCluster(int nPeaks); // Not Used for now


 public:
  BeamGEMProjection();
  BeamGEMProjection(TString projName, Int_t nch);
  ~BeamGEMProjection();
  
  inline vector< AHit> GetHits() const {return vHits;};
  inline Double_t GetChargeSum() const { return charge_sum;};
  
  inline int GetNHits() const {return nHits;};
  inline TString GetProjName() const{return strProjName;};
  inline int GetNStrips() const {return nStrips;};
  inline TH1D* GetRawHist() const {return h_raw;};
  inline TH1D* GetProjHist() const {return h_proj;};
  
  inline double GetBaselineMean() const {return baseline_mean;};
  inline double GetBaselineRMS() const {return baseline_rms;};
  inline double GetOverallMean() const {return overall_mean;};
  inline double GetOverallRMS() const {return overall_rms;};

  // Called by Users
  int Process();

  void AddStrip(BeamGEMStrip* );
  void PlotResults(TString, int);
    
  //called by BeamGEMPlane::Process(), during post-check
  int PostProcess();
  void UpdateHits( vector<int> vHitsMask );

  ClassDef(BeamGEMProjection,0);
};
