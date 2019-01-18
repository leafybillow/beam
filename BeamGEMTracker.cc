#include "BeamGEMTracker.h"
#include "BeamGEMPlane.h"
#include "BeamGEMProjection.h"
#include "TH1D.h"
ClassImp(BeamGEMTracker);

BeamGEMTracker::BeamGEMTracker()
  :fSlope_zx(NULL),fSlope_zy(NULL),
   fTheta(NULL),fPhi(NULL),
   fDet_x(NULL),fDet_y(NULL),
   vPlanes(),
   nTracks(-1){
  
  Init();
}

BeamGEMTracker::~BeamGEMTracker(){

}
void BeamGEMTracker::Init(){


}


void BeamGEMTracker::PlotResults(TString runName, int ievt){
  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  TPad* pad_gem1_y = new TPad("pad_gem1_y","",0.0,1 ,0.6,0.5);
  TPad* pad_gem1_x = new TPad("pad_gem1_x","",0.6,1 ,1.0,0.5);
  TPad* pad_gem2_y = new TPad("pad_gem2_y","",0.0,0.5 ,0.6,0.0);
  TPad* pad_gem2_x = new TPad("pad_gem2_x","",0.6,0.5 ,1.0,0.0);
  c1->cd();
  pad_gem1_y->Draw();
  pad_gem1_x->Draw();
  pad_gem2_y->Draw();
  pad_gem2_x->Draw();
  
  pad_gem1_y->cd();
  vPlanes[0]->GetProjectionY()->GetTH1D()->Draw();
  pad_gem1_x->cd();
  vPlanes[0]->GetProjectionX()->GetTH1D()->Draw();

  pad_gem2_y->cd();
  vPlanes[1]->GetProjectionY()->GetTH1D()->Draw();
  pad_gem2_x->cd();
  vPlanes[1]->GetProjectionX()->GetTH1D()->Draw();

  c1->SaveAs(Form("%s-Tracker-evt-%d.png",
		  runName.Data(),
		  ievt));
  delete c1;
}

void BeamGEMTracker::AddPlane(BeamGEMPlane* bgPlane){
  vPlanes.push_back(bgPlane);
}
