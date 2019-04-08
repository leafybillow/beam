#include "BeamGEMData.h"

ClassImp(BeamGEMData);

BeamGEMData::BeamGEMData()
  :adc(),
   adc_sum(),
   common_mode(),
   id_strip(),
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


    for(int isample=0;isample<6;isample++){
      adc[isample] = new Double_t[ size ];
    }

    nChannel = -1; // nChannel will be overwritten when reading a branch
}

BeamGEMData::~BeamGEMData(){
  // Nothing Here
}

Int_t BeamGEMData::FindPeakTime(){
  Double_t channel_sum[6];
  for(int isample =0; isample<6;isample++){
    channel_sum[isample] = 0;
    for(int ich=0;ich<nChannel;ich++){
      channel_sum[isample]+=adc[isample][ich];
    }
  }

  Double_t peak_val = channel_sum[0];
  Double_t peak_sample = 0;
  int isample = 1;
  while( isample<6 ){

    if(channel_sum[isample]>peak_val){
      peak_val = channel_sum[isample];
      peak_sample = isample;
    }
    isample++;
  }
  return peak_sample;
}
