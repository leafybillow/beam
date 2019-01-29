#include "BeamGEMTracker.h"
#include "BeamGEMPlane.h"
#include "BeamGEMProjection.h"
#include "TH1D.h"
#include "TText.h"
ClassImp(BeamGEMTracker);

BeamGEMTracker::BeamGEMTracker()
  :fSlope_zx(),fSlope_zy(),
   fTheta(),fPhi(),
   fDet_x(),fDet_y(),
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
  TPad* pad_gem1_y = new TPad("pad_gem1_y","",0.0,0.9 ,0.6,0.45);
  TPad* pad_gem1_x = new TPad("pad_gem1_x","",0.6,0.9 ,1.0,0.45);
  TPad* pad_gem2_y = new TPad("pad_gem2_y","",0.0,0.45 ,0.6,0.0);
  TPad* pad_gem2_x = new TPad("pad_gem2_x","",0.6,0.45 ,1.0,0.0);
  c1->cd();
  pad_gem1_y->Draw();
  pad_gem1_x->Draw();
  pad_gem2_y->Draw();
  pad_gem2_x->Draw();
  
  pad_gem1_y->cd();
  vPlanes[0]->GetProjectionY()->GetHistogram()->Draw();
  pad_gem1_x->cd();
  vPlanes[0]->GetProjectionX()->GetHistogram()->Draw();

  pad_gem2_y->cd();
  vPlanes[1]->GetProjectionY()->GetHistogram()->Draw();
  pad_gem2_x->cd();
  vPlanes[1]->GetProjectionX()->GetHistogram()->Draw();

  c1->cd();
  TText *text= new TText(0.0,0.95,
			 Form("%s-Tracker-evt-%d",runName.Data(),ievt));
  text->Draw("same");
  
  c1->SaveAs(Form("%s-Tracker-evt-%d.png",
		  runName.Data(),
		  ievt));
  delete c1;
}

void BeamGEMTracker::AddPlane(BeamGEMPlane* bgPlane){
  vPlanes.push_back(bgPlane);
}
