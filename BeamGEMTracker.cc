#include "BeamGEMTracker.h"
#include "BeamGEMPlane.h"
#include "BeamGEMProjection.h"
#include "BeamParameters.h"
#include "BeamConfig.h"

#include "TH1D.h"
#include "TText.h"
#include "TBox.h"
#include "TGraph.h"
#include "TLinearFitter.h"
#include "TF1.h"
#include "TMath.h"

ClassImp(BeamGEMTracker);

BeamGEMTracker::BeamGEMTracker()
  :fSlope_zx(),fSlope_zy(),
   fTheta(),fPhi(),
   fDet_x(),fDet_y(),fDet_z(),
   fGEM_z(),fHit_x(),fHit_y(),
   qdc_value(0),
   vPlanes(),vTracks(),
   nTracks(0),isGoldenTrack(0),track_npt(0)
{
  lf = new TLinearFitter();
  lf->SetFormula("pol1");
}

BeamGEMTracker::~BeamGEMTracker(){
}

void BeamGEMTracker::Init(){
  nPlanes = vPlanes.size();
  
  vector <double> vec_buff;
  effNhits.clear();
  
  for(int i=0;i<nPlanes;i++){
    Int_t nhitsx = vPlanes[i]->GetProjectionX()->GetNHits() ;
    Int_t nhitsy = vPlanes[i]->GetProjectionY()->GetNHits() ;

    Int_t eff = (nhitsy>=nhitsx ? nhitsx:nhitsy);
    
    if(eff> 0 ){
      track_npt ++;
      effNhits.push_back(eff);      
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
  if(track_npt<=1)
    nTracks = 0;
  else if(track_npt==2){
    isGoldenTrack = 0;
    for(int i=0; i< effNhits[0]; i++){
      ATrack goodTrack;
      double min_slope = 100000; // some non-sense large number
      for(int j=0; j < effNhits[1];j++){
	int pattern[2] = {i,j};
	ATrack aTrack = GenerateCandidates(pattern);
	FitATrack(&aTrack);
	double slope = aTrack.fSlope_zy * aTrack.fSlope_zx;
	if(slope<min_slope){
	  min_slope = slope;
	  goodTrack =aTrack;
	}
      }
      if(effNhits[1]!=0)
	vTracks.push_back(goodTrack);	  
      if(effNhits[1]>1){
	int myid = goodTrack.myPattern[1];
	SwapHits(1, myid, effNhits[1]-1);
	effNhits[1] = effNhits[1] -1;
      }
      else if(effNhits[1] ==1)
	effNhits[1] = 0;
    }
  }
  else if(track_npt>2){
    isGoldenTrack = 1;
    int imax = effNhits[0];
    int jmax = effNhits[1];
    for(int i=0; i<imax ; i++){
      vector<ATrack> vTracks_holder;
      for(int j=0; j < jmax;j++){
	
	ATrack this_track = PingForward(i,j);
	
	if(this_track.fChi2>0)
	  vTracks_holder.push_back(this_track);
	
      }// End of second plane loop
      
      if(vTracks_holder.size()!=0){
	vector<ATrack>::iterator itk = vTracks_holder.begin();
	vector<ATrack>::iterator itk_next = itk+1;
	while(itk_next!=vTracks_holder.end()){
	  if( (*itk_next).fChi2< (*itk).fChi2 )
	    itk  = itk_next;
	  itk_next++;
	}
	
	vTracks.push_back(*itk);
	//swap out picked points
	int myID =  ((*itk).myPattern)[1];
	SwapHits(1,myID,jmax-1);
	jmax= jmax-1;
	if(effNhits[2]>1){
	  myID = ((*itk).myPattern)[2];
	  SwapHits(2,myID,effNhits[2]-1);
	  effNhits[2] = effNhits[2] -1;
	}
	else if(effNhits[2]==1)
	  effNhits[2] = effNhits[2] -1;
      }
    } // End of first plane loop
  } // End of else if track_np>2 
  
  vector<ATrack>::iterator it = vTracks.begin();
  while(it!=vTracks.end()){
    FitATrack(&(*it));
    nTracks++;
    it++;
  }
  ProjectHits(); // To detector Planes
}

void BeamGEMTracker::ProjectHits(){
  vector<double>::iterator idetz = fDet_z.begin();
  double slope_zx, slope_zy;
  double intercept;
  
  while(idetz!=fDet_z.end()){
    vector<double> det_hitx;
    vector<double> det_hity;
    vector<double> det_theta;
    vector<double> det_phi;
    for(int i=0;i<nTracks;i++){

      double z = *idetz;

      slope_zx = vTracks[i].fSlope_zx;
      intercept = vTracks[i].fIntercept_x;
      double x = slope_zx*z +intercept;
      slope_zy = vTracks[i].fSlope_zy;
      intercept = vTracks[i].fIntercept_y;
      double y = slope_zy*z +intercept;
      det_hity.push_back(y);
      det_hitx.push_back(x);

      double tan_theta = TMath::Sqrt( slope_zx*slope_zx + slope_zy* slope_zy);
      double theta = TMath::ATan( tan_theta) * TMath::RadToDeg();
      if (slope_zy>0)
	theta = -theta;
      det_theta.push_back(theta);

      double tan_phi = slope_zy/slope_zx;
      double phi = TMath::ATan(tan_phi)*TMath::RadToDeg();
      det_phi.push_back(phi);
    
    }
    fDet_x.push_back(det_hitx);
    fDet_y.push_back(det_hity);
    
    fDet_theta.push_back(det_theta);
    fDet_phi.push_back(det_phi);
    idetz++;
  }
}

ATrack BeamGEMTracker::PingForward(int i, int j){ // hits id

  // Project from i to j
  
  double distance_cut = 10; // (mm)
  double distance = distance_cut; // an non-sense large number
  
  double z1 = fGEM_z[0];
  double z2 = fGEM_z[1];
  double z3 = fGEM_z[2];

  double x1 = fHit_x[0][i];
  double x2 = fHit_x[1][j];
  double x3 = (x1-x2)/(z1-z2)*(z3-z1)+x1;
    
  double y1 = fHit_y[0][i];
  double y2 = fHit_y[1][j];
  double y3 = (y1-y2)/(z1-z2)*(z3-z1)+y1;

  ATrack aTrack;
  if(fabs(x3)>50 || fabs(y3)>100 ){ // out of boundary
    int pattern[3] = {i,j,0}; // the last index is a dummy
    aTrack = GenerateCandidates(pattern);
    aTrack.fChi2 = -1; // an invalid track
    return aTrack;
  }
  int nhits = effNhits[2];
  int idFound=0;
  int isFound=0;
  for(int ihit=0; ihit<nhits;ihit++){
    
    double delta_x = fabs(fHit_x[2][ihit] - x3);
    double delta_y = fabs(fHit_y[2][ihit] - y3);
    if(delta_x+delta_y<distance){
      idFound = ihit;
      distance = delta_x + delta_y;
      isFound = 1;
    }
  }
  if(isFound){
    int pattern[3] = {i,j,idFound};
    aTrack = GenerateCandidates(pattern);
    aTrack.fChi2 = distance;
  }
  else{
    int pattern[3] = {i,j,0};
    aTrack = GenerateCandidates(pattern);
    aTrack.fChi2 = -1;
  }
    
  return aTrack;
}
void BeamGEMTracker::SwapHits(int iplane, int i, int j){

  double buff_hitx = fHit_x[iplane][i];
  double buff_hity = fHit_y[iplane][i];

  fHit_x[iplane][i] = fHit_x[iplane][j];
  fHit_y[iplane][i] = fHit_y[iplane][j];

  fHit_x[iplane][j] = buff_hitx;
  fHit_y[iplane][j] = buff_hity;

}

ATrack BeamGEMTracker::GenerateCandidates(int* pattern){
  ATrack aTrack ;
  Int_t npt = fGEM_z.size();
  
  for(int ipt =0; ipt<npt; ipt++){
    (aTrack.z).push_back(fGEM_z[ipt]);
    (aTrack.x).push_back(fHit_x[ipt][ pattern[ipt] ]);
    (aTrack.y).push_back(fHit_y[ipt][ pattern[ipt] ]);

    (aTrack.myPattern).push_back( pattern[ipt] );
  }
  return aTrack;
}

bool BeamGEMTracker::FitATrack(ATrack* aTrack){
  
  Int_t npt = (aTrack->z).size();
  Double_t *z = new Double_t[npt]; // in z direction
  Double_t *x = new Double_t[npt];
  Double_t *y = new Double_t[npt];
  
    
  for(int i=0;i<npt;i++){
    z[i] = (aTrack->z)[i];
    y[i] = (aTrack->y)[i];
    x[i] = (aTrack->x)[i];
  }
  lf->ClearPoints();
  lf->AssignData(npt,1,z,y);
  lf->Eval();

  aTrack->fIntercept_y = lf->GetParameter(0);
  aTrack->fSlope_zy = lf->GetParameter(1);
  aTrack->fChi2 = lf->GetChisquare();
  
  lf->ClearPoints();
  lf->AssignData(npt,1,z,x);
  lf->Eval();

  aTrack->fIntercept_x = lf->GetParameter(0);
  aTrack->fSlope_zx = lf->GetParameter(1);
  aTrack->fChi2 += lf->GetChisquare();

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

  vector< TBox* > box_det_zy;
  vector< TBox* > box_det_zx;
  
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

  int nDets = fDet_z.size();
  for(int iDet =0 ; iDet<nDets; iDet++){
    double pos_z = fDet_z[iDet];
    double thickness = fDet_thickness[iDet];
    double width_x = fDet_width_x[iDet];
    double width_y = fDet_width_y[iDet];

    double pos_x = fDet_pos_x[iDet];
    double pos_y = fDet_pos_y[iDet];

    TBox *bDet_zy = new TBox( pos_z-thickness/2.0, pos_y - width_y/2.0,
			      pos_z+thickness/2.0, pos_y + width_y/2.0);
    TBox *bDet_zx = new TBox( pos_z-thickness/2.0, pos_x - width_x/2.0,
			      pos_z+thickness/2.0, pos_x + width_x/2.0);
    bDet_zx->SetFillColor(kRed);
    bDet_zx->SetFillStyle(3144);
    bDet_zy->SetFillColor(kRed);
    bDet_zy->SetFillStyle(3144);


    box_det_zx.push_back(bDet_zx);
    box_det_zy.push_back(bDet_zy);
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

    for(int idet=0; idet<nDets;idet++){
      pad_track->cd(1);
      box_det_zy[idet]->Draw("same");
      pad_track->cd(2);
      box_det_zx[idet]->Draw("same");

    }
  }
  
  TF1 *flin_zx;
  TF1 *flin_zy;
  
  vector<TF1*> vlin_zx;
  vector<TF1*> vlin_zy;
  
  double par[2];  
  for(int i=0;i<nTracks;i++){
    flin_zx = new TF1(Form("fzx%d",i),"pol1",-10,10e3);
    flin_zy = new TF1(Form("fzy%d",i),"pol1",-10,10e3);

    par[1] = vTracks[i].fSlope_zy;
    par[0] = vTracks[i].fIntercept_y;
    flin_zy->SetParameters(par);

    if(track_npt == 2)
      flin_zy->SetLineStyle(2);
    vlin_zy.push_back(flin_zy);
    
    par[1] = vTracks[i].fSlope_zx;
    par[0] = vTracks[i].fIntercept_x;
    flin_zx->SetParameters(par);

    if(track_npt == 2)
      flin_zx->SetLineStyle(2);
    vlin_zx.push_back(flin_zx);
  }
  for(int i=0;i<nTracks;i++){
    pad_track->cd(1);
    vlin_zy[i]->Draw("same");
    pad_track->cd(2);
    vlin_zx[i]->Draw("same");
  }
  
  c1->cd();
  TString qdc_text = "QDC:";
  qdc_text += Form("%d \t ", (int)qdc_value);
  
  TText *text= new TText(0.0,0.95,
			 Form("%s-Tracker-evt-%d, %s",
			      runName.Data(),ievt, qdc_text.Data() ));
  text->Draw("same");
  
  c1->SaveAs(Form("%s-Tracker-evt-%d.png",
		  runName.Data(),
		  ievt));
  delete c1;
}

void BeamGEMTracker::LoadDetectorGeometry(BeamConfig *fConfig){
 
  fDet_z = fConfig->GetZ_Det();
  fDet_pos_x = fConfig->GetPositionX_Det();
  fDet_pos_y = fConfig->GetPositionY_Det();

  fDet_width_x = fConfig->GetWidthX_Det();
  fDet_width_y = fConfig->GetWidthY_Det();
  fDet_thickness = fConfig->GetThickness_Det();

}
