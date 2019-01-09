#include <TObject.h>
#include <RooInt.h>
#include <vector>
#include <TH1.h>
using namespace std;

class BeamGEMPlane: public TObject{
 private:
  vector<double> fPos_x; //position in GEM coordinate, unit: mm
  vector<double> fPos_y;
  vector<double> fCharge_x; // charge amplitude, unit: adc
  vector<double> fCharge_y;
  vector<double> fCorelation; // (x-y)/(x+y)

  int nHits;  //Number of hits found

 public:
  BeamGEMPlane();
  ~BeamGEMPlane();

  /* inline void SetPositionX(vector<double> pos){pos_x = pos;}; */
  /* inline void SetPositionY(vector<double> pos){pos_y = pos;}; */
  /* inline void SetChargeX(vector<double> charge){charge_x = charge;}; */
  /* inline void SetChargeY(vector<double> charge){charge_y = charge;}; */

  /* inline void  SetNHits(int n){nhits= n; }; */


  inline vector<double> GetPositionX() const {return fPos_x;};
  inline vector<double> GetPositionY() const {return fPos_y;};
  inline int GetNHits() const {return nHits;};

  void CheckHits();  
  void Init();
  void Process();
  ClassDef(BeamGEMPlane,0);
};
