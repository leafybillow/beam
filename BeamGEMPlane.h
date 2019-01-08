#include <TObject.h>
#include <RooInt.h>
#include <vector>
#include <TH1.h>
using namespace std;

class BeamGEMPlane: public TObject{
 private:
  vector<double> pos_x; //position in GEM coordinate, unit: mm
  vector<double> pos_y;
  vector<double> charge_x; // charge amplitude, unit: adc
  vector<double> charge_y;
  vector<double> res_x; // spatial resolution, unit: um
  vector<double> res_y;
  vector<double> corelation; // (x-y)/(x+y)
  vector<int> mpl_x; //multiplicity 
  vector<int> mpl_y;
  
  TH1D *h1d_x; // one dimensional historgram for charge distribution
  TH1D *h1d_y;

  double cmn[6]; // common mode noise on each APC
  int nhits;  //Number of hits found
  bool isSplit; // Does it has splitting peaks ?

 public:
  BeamGEMPlane();
  ~BeamGEMPlane();

  inline void SetPositionX(vector<double> pos){pos_x = pos;};
  inline void SetPositionY(vector<double> pos){pos_y = pos;};
  inline void SetChargeX(vector<double> charge){charge_x = charge;};
  inline void SetChargeY(vector<double> charge){charge_y = charge;};

  inline void  SetNHits(int n){nhits= n; };
  
  inline vector<double> GetPositionX() const {return pos_x;};
  inline vector<double> GetPositionY() const {return pos_y;};
  inline int GetNHits(){return nhits;};
  
  ClassDef(BeamGEMPlane,0);
};
