#include <TObject.h>
#include <RooInt.h>
#include <vector>
using namespace std;

class BeamGEMTracker: public TObject{
 private:
  vector<double> pos_x;
  vector<double> pos_y;
  int nhits;
 public:
  BeamGEMTracker();
  ~BeamGEMTracker();

  void SetPositionX(vector<double> pos);
  void SetPositionY(vector<double> pos);
  void SetNHits(int n);
  vector<double> GetPositionX();
  vector<double> GetPositionY();
  int GetNHits();
  
  ClassDef(BeamGEMTracker,0);
};
