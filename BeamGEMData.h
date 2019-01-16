#include <TObject.h>
#include <RooInt.h>

class BeamGEMData : public TObject{
  
 private:
 public:
  BeamGEMData();
  BeamGEMData(Int_t proj_id, Int_t size);
  ~BeamGEMData();
  
  // containers for data
  Double_t* adc[6];
  Double_t* adc_sum;
  Double_t* common_mode;
  Double_t* id_strip; // channel-to-strip mapping array
  Double_t nChannel;  // Note: nch in SBS-offline decoder is a float number

  // identification for projection
  Int_t id_Proj; // starts from 0 to 3; {"x1","y1","x2","y2"}

  ClassDef(BeamGEMData,0);
};
