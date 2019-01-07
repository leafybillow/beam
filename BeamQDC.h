#include <TObject.h>
#include <RooInt.h>

class BeamQDC: public TObject{
  
 public: 

  //Essentially treat it as a C structure
  double us_lo; // upstream detector low range data
  double us_hi; // high range data
  double ds_lo; // Dowstream detector
  double ds_hi;

  BeamQDC();
  ~BeamQDC();

  ClassDef(BeamQDC,0);
};
