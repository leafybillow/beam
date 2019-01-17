#include "BeamGEMPlane.h"
#include <iostream>
#include "BeamGEMProjection.h"

ClassImp(BeamGEMPlane);

BeamGEMPlane::BeamGEMPlane(TString name)
  :fPos_x(NULL),fPos_y(NULL),
   fCharge_x(NULL),fCharge_y(NULL),
   fWidth_x(NULL),fWidth_y(NULL),
   fCorelation(NULL),
   vHitsMask_x(NULL),vHitsMask_y(NULL),
   nHits(-1),
   bgProjX(NULL),bgProjY(NULL)
{  
  strPlaneName = name;
}
BeamGEMPlane::~BeamGEMPlane(){}


int BeamGEMPlane::Process(){
  int status = CheckProjections();
  if(status==0){

    nHits = Reconstruct();
    if(nHits == 1){
      bgProjY->UpdateHits(vHitsMask_y);
      bgProjX->UpdateHits(vHitsMask_x);
      CollectResults();
    }

  }
  else
    return 1; // Fail PlaneName Check
}

int BeamGEMPlane::CheckHits(){
  // Final Check 
  int nHitsX = fPos_x.size();
  int nHitsY = fPos_y.size();
  if(nHitsY!=nHitsX){
    std::cerr << "Error: "
	      << __FILE__ << ":"
	      << __LINE__ << ":"
	      << __FUNCTION__  << ":"
	      << "Number of Hits between Projection mismatch"
	      <<std::endl;

    return 1;
  }
  else
    return 0;
}

int BeamGEMPlane::CheckProjections(){
  TString strName_Y = bgProjY->GetProjName();
  TString strName_X = bgProjX->GetProjName();
  if(strName_Y.Contains("y") && strName_X.Contains("x"))
    return 0;
  else{
    std::cerr << "Error:"
	      <<__FUNCTION__
	      <<" Incorrect Projections pair"
	      << std::endl;
    return 1;
  }
}

int BeamGEMPlane::Reconstruct(){
  return Reconstruct_v0();
}

double BeamGEMPlane::EvalCorrelation(double x, double y){
  //Experimenting Function
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
  
  c1->SaveAs(Form("%s-%s-evt-%d.pdf",
		  runName.Data(),
		  strPlaneName.Data(),
		  ievt));
  
  delete c1;
}

void BeamGEMPlane::AddProjectionX(BeamGEMProjection* bgProj){
  bgProjX = bgProj;
}

void BeamGEMPlane::AddProjectionY(BeamGEMProjection* bgProj){
  bgProjY = bgProj;
}

void BeamGEMPlane::CollectResults(){
  vector< AHit> vHits_x = bgProjX->GetHits();
  vector< AHit> vHits_y = bgProjY->GetHits();

  int nHits_x = vHits_x.size();
  int nHits_y = vHits_y.size();

  for(int iHits=0;iHits<nHits_x;iHits++){
    fCharge_x.push_back( vHits_x[iHits].fCharge);
    fPos_x.push_back( vHits_x[iHits].fPosition);
    fWidth_x.push_back( vHits_x[iHits].fWidth);
  }

  for(int iHits=0;iHits<nHits_y;iHits++){
    fCharge_y.push_back( vHits_y[iHits].fCharge);
    fPos_y.push_back( vHits_y[iHits].fPosition);
    fWidth_y.push_back( vHits_y[iHits].fWidth);
  }
}

int BeamGEMPlane::Reconstruct_v0(){
  // ver 0
  // Only used to separate out zero and single hit

  int nHits_x = bgProjX->GetNHits();
  int nHits_y = bgProjY->GetNHits();

  vHitsMask_x.clear();
  vHitsMask_y.clear();
  if(nHits_y==0 ||nHits_x==0){
    return 0;  
  }
  else if(nHits_y==1 && nHits_x==1){
    vHitsMask_x.push_back(1);
    vHitsMask_y.push_back(1); 
    return 1;      // accept only one hit
  }
  else if( nHits_y ==1 && nHits_x ==2) {
    vHitsMask_y.push_back(1); 
    vHitsMask_x.push_back(1);
    if( bgProjX->TestCrossTalk(0,1) ){
      vHitsMask_x.push_back(0);
      return 1 ; 
    }
    else
      return -1 ; // Not ready to resolve this type of hits yet
  }
  else if(nHits_x ==1 && nHits_y ==2 ){  
    vHitsMask_y.push_back(1); 
    vHitsMask_x.push_back(1);
    if( bgProjY->TestCrossTalk(0,1) ){
      vHitsMask_y.push_back(0);
      return 1 ; 
    }
    else
      return -1 ;
  }
  else if ( nHits_y ==2 && nHits_x == 2){
    vHitsMask_y.push_back(1); 
    vHitsMask_x.push_back(1);
    if( bgProjY->TestCrossTalk(0,1) &&
	bgProjX->TestCrossTalk(0,1) ){

      vHitsMask_y.push_back(0);
      vHitsMask_x.push_back(0);
      return 1 ; 
    }
    else
      return -1 ;
  }
  else 
    return -1;  // ignore this events
}

