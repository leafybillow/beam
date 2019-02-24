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
   nTracks(-1),isGoldenTrack(0){
  
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
  TH1D *h_buff;
  
  pad_gem1_y->cd();
  h_buff =vPlanes[0]->GetProjectionY()->GetRawHist();
  h_buff->Draw();
  h_buff =vPlanes[0]->GetProjectionY()->GetProjHist();
  h_buff->SetFillColor(kGreen);
  h_buff->SetLineColor(kWhite);
  h_buff->SetLineWidth(0.0);
  h_buff->Draw("box same");
  h_buff =vPlanes[0]->GetProjectionY()->GetRawHist();
  h_buff->Draw("same");  

  pad_gem1_x->cd();
  h_buff =vPlanes[0]->GetProjectionX()->GetRawHist();
  h_buff->Draw();
  h_buff =vPlanes[0]->GetProjectionX()->GetProjHist();
  h_buff->SetFillColor(kGreen);
  h_buff->SetLineColor(kWhite);
  h_buff->SetLineWidth(0.0);
  h_buff->Draw("box same");
  h_buff =vPlanes[0]->GetProjectionX()->GetRawHist();
  h_buff->Draw("same");

  pad_gem2_y->cd();
  h_buff =vPlanes[1]->GetProjectionY()->GetRawHist();
  h_buff->Draw();
  h_buff =vPlanes[1]->GetProjectionY()->GetProjHist();
  h_buff->SetFillColor(kGreen);
  h_buff->SetLineColor(kWhite);
  h_buff->Draw("box same");
  h_buff =vPlanes[1]->GetProjectionY()->GetRawHist();
  h_buff->Draw("same");

  pad_gem2_x->cd();
  h_buff =vPlanes[1]->GetProjectionX()->GetRawHist();
  h_buff->Draw();
  h_buff =vPlanes[1]->GetProjectionX()->GetProjHist();
  h_buff->SetFillColor(kGreen);
  h_buff->SetLineColor(kWhite);
  h_buff->Draw("box same");
  h_buff =vPlanes[1]->GetProjectionX()->GetRawHist();
  h_buff->Draw("same");

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

void BeamGEMTracker::Process(){
  // Int_t nplane = vPlanes.size();
  // if(nplane==2){

  //   Int_t nhits1 = vPlanes[0]->GetNHits();
  //   Int_t nhits2 = vPlanes[1]->GetNHits();
  //   if(nhits2==1 && nhits1==1){
      

  //   }
    
}


// }

// void BeamGEMTracker::Tracking_v2(){
  

// }

// void BeamGEMTracker::Tracking_v3(){

  
// }
