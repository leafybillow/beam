#include "BeamGEMStrip.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

ClassImp(BeamGEMStrip);

BeamGEMStrip::BeamGEMStrip(double* d){
  fAmpl_raw = 0.0;
  fAmpl_fit = 0.0;
  fTau = 0.0;
  fT_start =0.0;
  fT_min = 0.0;

  WriteSamples(d);
  Init();
}

BeamGEMStrip::~BeamGEMStrip(){

}

void BeamGEMStrip::WriteSamples(double* d){
  // Number of sample is hardcoded here.
  // Probably we need to check size of array before Write data
  for(int i=0;i<6;i++){
    fData[i] = d[i];
  }
}

double BeamGEMStrip::ADCSum(){
  double ret =0;
  for(int i=0;i<6;i++){
    ret += fData[i];
  }
  return ret;
}

int BeamGEMStrip::FindMinimum(){
  double min = fData[0]; 
  int t_min = 0;
  for(int i=1;i<6;i++){
    if(fData[i]<min){
      min = fData[i];
      t_min = i;
    }
  }
  return t_min;
}

void BeamGEMStrip::FitData(){
  // Fit the histogram for now. May want to use chi2 fit without calling a histogram
  TH1D *h_fit =  new TH1D("h_fit","histogram for fit",6,-0.5,5.5);
  for(int i=0;i<6;i++){
    h_fit->SetBinContent(i+1,fData[i]);
  }
  TF1* fcn_sig = new TF1("fcn_sig",CRRCShaping,-10.0,10.0,4);

  double par[4];
  par[2] = 2.0; // An initial Guess, it means 2*25ns
  par[1]= fT_min-par[2];
  par[3]= 0.0; // After pedestal and common mode corrections ,a zero pedestal is expected.
  par[0]= fAmpl_raw*2.718;
  fcn_sig->SetParameters(par);
  
  fcn_sig->FixParameter(3,0.0); // Test: Fix pedestal to zero
  
  h_fit->Fit("QN","",0,6);

  fcn_sig->GetParameters(par);

  fAmpl_fit = par[0];
  fTau = par[2];
  fT_start = par[1];
}

void BeamGEMStrip::Process(){
  FitData();
  //PlotFitResult();
}

void BeamGEMStrip::Init(){
  fT_min = FindMinimum();
  fAmpl_raw = fData[fT_min];
  fIntegral = ADCSum();
}

double BeamGEMStrip::CRRCShaping(double* x, double* par){

  double t = x[0];
  double V_0 = par[0];
  double t_start = par[1];
  double tau=par[2];
  double offset=par[3];

  double fcn_val;

  if(t<t_start)
    fcn_val =offset;
  else
    fcn_val = -V_0*(t-t_start)/tau*TMath::Exp(-(t-t_start)/tau) + offset;

  return fcn_val;
}
