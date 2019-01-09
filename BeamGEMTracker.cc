#include "BeamGEMTracker.h"
#include "BeamGEMPlane.h"

ClassImp(BeamGEMTracker);

BeamGEMTracker::BeamGEMTracker(BeamGEMPlane* planeX, BeamGEMPlane* planeY)
  :fSlope_zx(NULL),fSlope_zy(NULL),fTheta(NULL),fPhi(NULL),fDet_x(NULL),fDet_y(NULL),nTracks(-1){
  
  Init();
}

BeamGEMTracker::~BeamGEMTracker(){

}
