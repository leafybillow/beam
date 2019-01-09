#include "BeamGEMProjection.h"
ClassImp(BeamGEMProjection);

BeamGEMProjection::BeamGEMProjection(TH1D* h, double *rms)
  :fPosition(NULL),fCharge(NULL),fRes(NULL),fMpl(NULL),nStrips(-1),nHits(-1),isSplit(0){
  
  h_proj =h;
  Init();
}

BeamGEMProjection::~BeamGEMProjection(){

}


void BeamGEMProjection::Init(){
  nStrips = h_proj->GetNbinsX();
}

void BeamGEMProjection::Process(){


}




