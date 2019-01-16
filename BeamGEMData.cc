#include "BeamGEMData.h"

ClassImp(BeamGEMData);

BeamGEMData::BeamGEMData()
  :adc(),
   adc_sum(NULL),
   common_mode(NULL),
   id_strip(NULL),
   nChannel(-1),
   id_Proj(-1)
{

}
BeamGEMData::BeamGEMData(Int_t proj_id, Int_t size)
{
    id_strip = new Double_t[ size ];
    adc_sum = new Double_t[ size ];
    common_mode = new Double_t[ size ];
    id_Proj = proj_id;


    for(int iadc=0;iadc<6;iadc++){
      adc[iadc] = new Double_t[ size ];
    }

    nChannel = -1; // nChannel will be overwritten when reading a branch
}

BeamGEMData::~BeamGEMData(){
  // Nothing Here
}
