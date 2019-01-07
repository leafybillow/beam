#include "BeamGEMTracker.h"
ClassImp(BeamGEMTracker);

BeamGEMTracker::BeamGEMTracker()
  :pos_x(NULL),nhits(0){

}

BeamGEMTracker::~BeamGEMTracker(){

}

void BeamGEMTracker::SetPositionX(vector<double> pos){
  pos_x = pos;
}

void BeamGEMTracker::SetNHits(int n){
  nhits = n;
}

vector<double> BeamGEMTracker::GetPositionX(){
  return pos_x;
}
