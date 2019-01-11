#include "BeamGEMPlane.h"
#include <iostream>
#include "BeamGEMProjection.h"

ClassImp(BeamGEMPlane);

BeamGEMPlane::BeamGEMPlane()
  :fPos_x(NULL),fPos_y(NULL),fCharge_x(NULL),fCharge_y(NULL),fCorelation(NULL),nHits(-1),bgProjX(NULL),bgProjY(NULL){

}
BeamGEMPlane::~BeamGEMPlane(){}

int BeamGEMPlane::Process(){
  int status = CheckProjections();
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

int BeamGEMPlane::CheckProjections(){
  TString strName_Y = bgProjY->GetProjName();
  TString strName_X = bgProjX->GetProjName();
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

void BeamGEMPlane::MatchHits(){
  
}

double BeamGEMPlane::EvalCorrelation(double x, double y){

  double ret;
  ret = (x-y)/(x+y);
  return ret;
}

void BeamGEMPlane::PlotResults(TString runName, int ievt){

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(1,2);

  c1->cd(1);
  bgProjY->GetTH1D()->Draw();

  TVirtualPad *c2 = c1->cd(2);
  c2->Divide(2,1);
  c2->cd(1);
  bgProjX->GetTH1D()->Draw();
  
  c1->SaveAs(Form("%s-evt%d-Plane.pdf",runName.Data(),ievt));
  
  delete c1;
}

void BeamGEMPlane::AddProjectionX(BeamGEMProjection* bgProj){
  bgProjX = bgProj;
}

void BeamGEMPlane::AddProjectionY(BeamGEMProjection* bgProj){
  bgProjY = bgProj;
}
