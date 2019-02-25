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
  TCanvas *c1 = new TCanvas("c1","c1",1400,1000);
  c1->Divide(2,1);
  TVirtualPad* pad_track = c1->cd(1);
  TVirtualPad* pad_gems = c1->cd(2);

  Int_t nplane =vPlanes.size();
  pad_gems->Divide(1,nplane);

  TH1D *h_buff;
  for(int iplane =0 ;iplane<nplane;iplane++){
    pad_gems->cd(iplane+1);
    TPad* pad_gem1_y = new TPad("pad_gem1_y","",0.0,1.0 ,0.66,0.0);
    TPad* pad_gem1_x = new TPad("pad_gem1_x","",0.66,1.0 ,1.0,0.0);
    pad_gem1_y->Draw();
    pad_gem1_x->Draw();
    
    pad_gem1_y->cd();
    h_buff =vPlanes[iplane]->GetProjectionY()->GetRawHist();
    h_buff->SetStats(0);
    h_buff->Draw();
    h_buff =vPlanes[iplane]->GetProjectionY()->GetProjHist();
    h_buff->SetStats(0);
    h_buff->SetFillColor(kGreen);
    h_buff->SetLineColor(kWhite);
    h_buff->SetLineWidth(0.0);
    h_buff->Draw("box same");
    h_buff =vPlanes[iplane]->GetProjectionY()->GetRawHist();
    h_buff->Draw("same");  

    pad_gem1_x->cd();
    h_buff =vPlanes[iplane]->GetProjectionX()->GetRawHist();
    h_buff->SetStats(0);
    h_buff->Draw();
    h_buff =vPlanes[iplane]->GetProjectionX()->GetProjHist();
    h_buff->SetStats(0);
    h_buff->SetFillColor(kGreen);
    h_buff->SetLineColor(kWhite);
    h_buff->SetLineWidth(0.0);
    h_buff->Draw("box same");
    h_buff =vPlanes[iplane]->GetProjectionX()->GetRawHist();
    h_buff->Draw("same");  

  }

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
