#include "BeamGEMTracker.h"
#include "BeamGEMPlane.h"
#include "BeamGEMProjection.h"
#include "TH1D.h"
#include "TText.h"
#include "TBox.h"
#include "TGraph.h"
#include "TLinearFitter.h"

ClassImp(BeamGEMTracker);

BeamGEMTracker::BeamGEMTracker()
  :fSlope_zx(),fSlope_zy(),
   fTheta(),fPhi(),
   fDet_x(),fDet_y(),fDet_z(),
   fGEM_z(),fHit_x(),fHit_y(),
   vPlanes(),
   nTracks(-1),isGoldenTrack(0),track_flag(1)
{
  lf = new TLinearFitter();
  lf->SetFormula("pol1");
}

BeamGEMTracker::~BeamGEMTracker(){
}

void BeamGEMTracker::Init(){
  nPlanes = vPlanes.size();
  
  vector <double> vec_buff;
  
  for(int i=0;i<nPlanes;i++){
    fGEM_z.push_back( vPlanes[i]->GetPositionZ() );
    track_flag *= vPlanes[i]->GetNHits() ;

    vec_buff = vPlanes[i]->GetPositionX();
    fHit_x.push_back(vec_buff);
    vec_buff = vPlanes[i]->GetPositionY();
    fHit_y.push_back(vec_buff);
  }
  
}

void BeamGEMTracker::Process(){
  Init();
  if(track_flag == 1)
    isGoldenTrack = FitSingleTrack(0);
  
}

bool BeamGEMTracker::FitSingleTrack(int iHit){
  Double_t *z = new Double_t[nPlanes]; // in z direction
  Double_t *x = new Double_t[nPlanes];
  Double_t *y = new Double_t[nPlanes];
  
  ATrack aTrack;
  
  for(int i=0;i<nPlanes;i++){
    z[i] = fGEM_z[i];
    y[i] = fHit_y[i][iHit];
    x[i] = fHit_x[i][iHit];
  }
  lf->AssignData(nPlanes,1,z,x);
  lf->Eval();

  aTrack.fIntercept = lf->GetParameter(0);
  aTrack.fSlope = lf->GetParameter(1);
  aTrack.fChi2 = lf->GetChisquare();
  vTrack_zx.push_back(aTrack);

  lf->AssignData(nPlanes,1,z,y);
  lf->Eval();

  aTrack.fIntercept = lf->GetParameter(0);
  aTrack.fSlope = lf->GetParameter(1);
  aTrack.fChi2 = lf->GetChisquare();
  vTrack_zy.push_back(aTrack);

  return 1;
}

void BeamGEMTracker::PlotResults(TString runName, int ievt){
  TCanvas *c1 = new TCanvas("c1","c1",1600,1000);
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
    h_buff->Draw("box same");
    h_buff =vPlanes[iplane]->GetProjectionX()->GetRawHist();
    h_buff->Draw("same");  

  }
  
  Double_t *gem_z = new Double_t[nplane];
  vector<Double_t> vec_hit_y;
  vector<Double_t> vec_hit_x;
  vector<Double_t> vec_pos_z;
  
  vector < TBox* > box_zy;
  vector < TBox* > box_zx;
  
  for(int i=0;i<nplane;i++){
    gem_z[i] = vPlanes[i]->GetPositionZ();
    int nhits_y = (vPlanes[i]->GetPositionY()).size();
    int nhits_x = (vPlanes[i]->GetPositionX()).size();
    if( nhits_x>0 && nhits_y >0){
      vec_hit_y.push_back( (vPlanes[i]->GetPositionY() )[0] );
      vec_hit_x.push_back( (vPlanes[i]->GetPositionX() )[0] );
      vec_pos_z.push_back( gem_z[i]);
    }
    
    TBox *bgem_y = new TBox(gem_z[i]+2,-100.0, gem_z[i]+20, 100);
    TBox *bgem_x = new TBox(gem_z[i]+2,-50.0, gem_z[i]+20, 50);
    
    bgem_y->SetLineColor(8);
    bgem_x->SetLineColor(8);
    box_zy.push_back(bgem_y);
    box_zx.push_back(bgem_x);
  }
  
  pad_track->Divide(1,2);
  pad_track->cd(1);
  int npt = vec_pos_z.size();
  Double_t *hit_y = new Double_t[npt];
  Double_t *hit_x = new Double_t[npt];
  Double_t *pos_z = new Double_t[npt];
  if(npt>0){
    for(int i=0;i<npt;i++){
      hit_y[i] = vec_hit_y[i];
      hit_x[i] = vec_hit_x[i];
      pos_z[i] = vec_pos_z[i];
    }
  
    TGraph* g_zy = new TGraph(npt, pos_z, hit_y);
    g_zy->SetMarkerSize(1);
    g_zy->SetMarkerStyle(47);
    g_zy->Draw("AP");
    g_zy->SetTitle("");
    g_zy->GetYaxis()->SetRangeUser(-110,110);
    g_zy->GetXaxis()->SetLimits(-10,850);
    
    for(int i=0;i<nplane;i++)
      box_zy[i]->Draw("l same");

  
    pad_track->cd(2);
    TGraph* g_zx = new TGraph(npt, pos_z, hit_x);
    g_zx->SetMarkerSize(1);
    g_zx->SetMarkerStyle(47);
    g_zx->Draw("AP");
    g_zx->SetTitle("");
    g_zx->GetYaxis()->SetRangeUser(-105,105);
    g_zx->GetXaxis()->SetLimits(-10,850);

    for(int i=0;i<nplane;i++)
      box_zx[i]->Draw("l same");
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

