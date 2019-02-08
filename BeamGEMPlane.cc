#include "BeamGEMPlane.h"
#include <iostream>
#include "BeamGEMProjection.h"

ClassImp(BeamGEMPlane);

BeamGEMPlane::BeamGEMPlane(TString name)
  :fPos_x(),fPos_y(),
   fCharge_x(),fCharge_y(),
   fWidth_x(),fWidth_y(),
   fSplit_x(),fSplit_y(),
   fCorelation(),
   vHitsMask_x(),vHitsMask_y(),
   nHits(-1),
   bgProjX(NULL),bgProjY(NULL)
{  
  strPlaneName = name;
}
BeamGEMPlane::~BeamGEMPlane(){}


int BeamGEMPlane::Process(){
  int status = CheckProjections();
  if(status==0){

    // nHits = Reconstruct();
    // if(nHits >0 ){
    //   bgProjY->UpdateHits(vHitsMask_y);
    //   bgProjX->UpdateHits(vHitsMask_x);
    //   CollectResults();
    // }
    int nHits_x = bgProjX->GetNHits();
    int nHits_y = bgProjY->GetNHits();
    nHits = (nHits_x >= nHits_y ? nHits_x : nHits_y);
    if(nHits == 1)
      CollectResults();
    return 0;
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
  
  return 0;
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
  bgProjY->GetHistogram()->Draw();

  TVirtualPad *c2 = c1->cd(2);
  c2->Divide(2,1);
  c2->cd(1);
  bgProjX->GetHistogram()->Draw();
  
  c1->SaveAs(Form("%s-%s-evt-%d.png",
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

  int nHits_x = bgProjX->GetNHits();
  int nHits_y = bgProjY->GetNHits();

  for(int iHits=0;iHits<nHits_x;iHits++){
    fCharge_x.push_back( vHits_x[iHits].fCharge);
    fPos_x.push_back( vHits_x[iHits].fPosition);
    fWidth_x.push_back( vHits_x[iHits].fWidth);
    // fSplit_x.push_back( vHits_x[iHits].fSplit);
  }

  for(int iHits=0;iHits<nHits_y;iHits++){
    fCharge_y.push_back( vHits_y[iHits].fCharge);
    fPos_y.push_back( vHits_y[iHits].fPosition);
    fWidth_y.push_back( vHits_y[iHits].fWidth);
    // fSplit_y.push_back( vHits_y[iHits].fSplit);
  }
}

