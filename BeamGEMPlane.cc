#include "BeamGEMPlane.h"
#include <iostream>
#include "BeamGEMProjection.h"

ClassImp(BeamGEMPlane);

BeamGEMPlane::BeamGEMPlane(BeamGEMProjection* projX, BeamGEMProjection* projY)
  :fPos_x(NULL),fPos_y(NULL),fCharge_x(NULL),fCharge_y(NULL),fCorelation(NULL),nHits(-1),gProjX(projX),gProjY(projY){

  Init();
}
BeamGEMPlane::~BeamGEMPlane(){}

void BeamGEMPlane::Init(){

}

int BeamGEMPlane::Process(){
  int status = CheckPlaneName();
  if(status==0){

    MatchHits();

    if(CheckHits()==0)
      return 0;
    else
      return 1; // Fail Number of hits Check
  }
  else
    return 1; // Fail PlaneName Check
}

int BeamGEMPlane::CheckHits(){

  int nHitsX = fPos_x.size();
  int nHitsY = fPos_y.size();
  if(nHitsY!=nHitsX){
    std::cout << "Error: "
	      <<__FUNCTION__ 
	      << " Number of Hits between Projection mismatch"
	      <<std::endl;

    return 1;
  }
  else
    return 0;
}

int BeamGEMPlane::CheckPlaneName(){
  TString strName_Y = gProjY->GetProjName();
  TString strName_X = gProjX->GetProjName();
  if(strName_Y=="y" && strName_X=="x")
    return 0;
  else{
    std::cout << "Error:"
	      <<__FUNCTION__
	      <<" Incorrect Projections pair"
	      << std::endl;
    return 1;
  }
}

void MatchHits(){
  
}

double EvalCorrelation(double x, double y){

  double ret;
  ret = (x-y)/(x+y);
  return ret;
}
