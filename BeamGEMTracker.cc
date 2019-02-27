#include "BeamGEMTracker.h"
#include "BeamGEMPlane.h"
#include "BeamGEMProjection.h"
#include "BeamParameters.h"

#include "TH1D.h"
#include "TText.h"
#include "TBox.h"
#include "TGraph.h"
#include "TLinearFitter.h"
#include "TF1.h"

ClassImp(BeamGEMTracker);

BeamGEMTracker::BeamGEMTracker()
  :fSlope_zx(),fSlope_zy(),
   fTheta(),fPhi(),
   fDet_x(),fDet_y(),fDet_z(),
   fGEM_z(),fHit_x(),fHit_y(),
   vPlanes(),vTraces(),
   nTracks(-1),isGoldenTrack(0),track_npt(0)
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
    Int_t nhitsx = vPlanes[i]->GetProjectionX()->GetNHits() ;
    Int_t nhitsy = vPlanes[i]->GetProjectionY()->GetNHits() ;
    
    if(nhitsx>0 && nhitsy>0){
      track_npt ++;
      
      fGEM_z.push_back( vPlanes[i]->GetPositionZ() );
      
      vec_buff = vPlanes[i]->GetPositionX();
      fHit_x.push_back(vec_buff);
      vec_buff = vPlanes[i]->GetPositionY();
      fHit_y.push_back(vec_buff);
    }
    
  }
  
}

void BeamGEMTracker::Process(){
  Init();
  GenerateCandidates();
  
  if(track_npt>=2){
    vector<ATrace>::iterator it = vTraces.begin();
    while(it!=vTraces.end()){
      FitATrack(&(*it));
      it++;
    }
    nTracks = vTraces.size();
    if(track_npt==3)
      isGoldenTrack=1;
  }
  else
    nTracks = 0;
}

void BeamGEMTracker::GenerateCandidates(){
  ATrace aTrace ;
  Int_t npt = fGEM_z.size();
  
  for(int ipt=0;ipt<npt;ipt++){
    (aTrace.z).push_back(fGEM_z[ipt]);
    (aTrace.x).push_back(fHit_x[ipt][0]);
    (aTrace.y).push_back(fHit_y[ipt][0]);
  }
  
  vTraces.push_back( aTrace );
  
}

bool BeamGEMTracker::FitATrack(ATrace* aTrace){
  
  Int_t npt = (aTrace->z).size();
  Double_t *z = new Double_t[npt]; // in z direction
  Double_t *x = new Double_t[npt];
  Double_t *y = new Double_t[npt];
  
    
  for(int i=0;i<npt;i++){
    z[i] = (aTrace->z)[i];
    y[i] = (aTrace->y)[i];
    x[i] = (aTrace->x)[i];
  }
  
  lf->AssignData(npt,1,z,y);
  lf->Eval();

  aTrace->fIntercept_y = lf->GetParameter(0);
  aTrace->fSlope_zy = lf->GetParameter(1);
  aTrace->fChi2 = lf->GetChisquare();
  
  lf->ClearPoints();
  lf->AssignData(npt,1,z,x);
  lf->Eval();

  aTrace->fIntercept_x = lf->GetParameter(0);
  aTrace->fSlope_zx = lf->GetParameter(1);
  aTrace->fChi2 += lf->GetChisquare();
  
  return 1;
}

void BeamGEMTracker::AddPlane(BeamGEMPlane* bgPlane){
  vPlanes.push_back(bgPlane);
}

void BeamGEMTracker::PlotResults(TString runName, int ievt){
  TCanvas *c1 = new TCanvas("c1","c1",1600,1000);
  c1->Divide(2,1);
  TVirtualPad* pad_track = c1->cd(1);
  TVirtualPad* pad_gems = c1->cd(2);

  pad_gems->Divide(1,nPlanes);

  TH1D *h_buff;
  for(int iplane =0 ;iplane<nPlanes;iplane++){
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
  
  vector<Double_t> vec_hit_y;
  vector<Double_t> vec_hit_x;
  vector<Double_t> vec_zpos_x;
  vector<Double_t> vec_zpos_y;
  
  vector < TBox* > box_zy;
  vector < TBox* > box_zx;

  for(int ipt=0;ipt<track_npt;ipt++){
    int nhitx = fHit_x[ipt].size();
    int nhity = fHit_y[ipt].size();

    for(int ihit=0; ihit<nhitx;ihit++){
      vec_hit_x.push_back( fHit_x[ipt][ihit] );
      vec_zpos_x.push_back( fGEM_z[ipt]);
    }

    for(int ihit=0; ihit<nhity;ihit++){
      vec_hit_y.push_back( fHit_y[ipt][ihit] );
      vec_zpos_y.push_back( fGEM_z[ipt]);
    }

  }
  
  for(int ipl=0;ipl<nPlanes;ipl++){
    double gem_z = vPlanes[ipl]->GetPositionZ();
    TBox *bgem_y = new TBox(gem_z+2,-100.0, gem_z+20, 100);
    TBox *bgem_x = new TBox(gem_z+2,-50.0, gem_z+20, 50);
    
    bgem_y->SetLineColor(8);
    bgem_x->SetLineColor(8);
    box_zy.push_back(bgem_y);
    box_zx.push_back(bgem_x);
  }
  
  pad_track->Divide(1,2);
  pad_track->cd(1);

  int nptx = vec_zpos_x.size();
  int npty = vec_zpos_y.size();

  Double_t *hit_x = new Double_t[nptx];
  Double_t *zpos_x = new Double_t[nptx];
  
  Double_t *hit_y = new Double_t[npty];
  Double_t *zpos_y = new Double_t[npty];
  // copy vector to array
  if(nptx>0 && npty>0){
    
    for(int i=0;i<nptx;i++){
      hit_x[i] = vec_hit_x[i];
      zpos_x[i] = vec_zpos_x[i];
    }

    for(int i=0;i<npty;i++){
      hit_y[i] = vec_hit_y[i];
      zpos_y[i] = vec_zpos_y[i];
    }

    TGraph* g_zy = new TGraph(npty, zpos_y, hit_y);
    g_zy->SetMarkerSize(2);
    g_zy->SetMarkerStyle(34);
    g_zy->Draw("AP");
    g_zy->SetTitle("");
    g_zy->GetYaxis()->SetRangeUser(-110,110);
    g_zy->GetXaxis()->SetLimits(-10,850);
    
    for(int i=0;i<nPlanes;i++)
      box_zy[i]->Draw("l same");

  
    pad_track->cd(2);
    TGraph* g_zx = new TGraph(nptx, zpos_x, hit_x);
    g_zx->SetMarkerSize(2);
    g_zx->SetMarkerStyle(34);
    g_zx->Draw("AP");
    g_zx->SetTitle("");
    g_zx->GetYaxis()->SetRangeUser(-105,105);
    g_zx->GetXaxis()->SetLimits(-10,850);

    for(int i=0;i<nPlanes;i++)
      box_zx[i]->Draw("l same");

  }
  
  TF1 *flin_zx;
  TF1 *flin_zy;
  double par[2];  
  for(int i=0;i<nTracks;i++){
    flin_zx = new TF1(Form("fzx%d",i),"pol1",-10,10e3);
    flin_zy = new TF1(Form("fzy%d",i),"pol1",-10,10e3);

    par[1] = vTraces[i].fSlope_zy;
    par[0] = vTraces[i].fIntercept_y;
    flin_zy->SetParameters(par);
    pad_track->cd(1);
    flin_zy->Draw("same");
    
    par[1] = vTraces[i].fSlope_zx;
    par[0] = vTraces[i].fIntercept_x;
    flin_zx->SetParameters(par);
    pad_track->cd(2);
    flin_zx->Draw("same");
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

