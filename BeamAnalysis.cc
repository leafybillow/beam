#include "BeamGEMData.h"
#include "BeamGEMStrip.h"
#include "BeamGEMPlane.h"
#include "BeamGEMProjection.h"
#include "BeamGEMTracker.h"

#include "BeamAnalysis.h"
#include "BeamConfig.h"

#include "TF1.h"
#include "TTree.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>

using namespace std;

#define N_GEM 2; // Number of GEMS

ClassImp(BeamAnalysis);

BeamAnalysis::BeamAnalysis(BeamConfig *beamConfig){
  fConfig = beamConfig;

  rf_raw = TFile::Open(fConfig->GetInputName());
  TString input_filename = rf_raw->GetName(); 
  cout << "--Input file " << input_filename << " is opened. " << endl;
  
  rf_output = TFile::Open(fConfig->GetOutputName() ,"RECREATE");
  TString output_filename = rf_output->GetName();
  cout << "--Output ROOTFile " << output_filename << " is recreated. " << endl;
}

BeamAnalysis::~BeamAnalysis(){
}

int BeamAnalysis::Process(){
  
  anaType = fConfig->GetAnalysisType();

  switch(anaType){
  case 0:
    Analysis(0);
    break;
  case 1:
    CalculatePed();
    break;
  case 2:
    CalculateRMS();
    break;
  case 3:
    Analysis(1);
    break;
  }
    
  rf_raw->Close();
  rf_output->Close();
  
  return 0;
}


int BeamAnalysis::CalculatePed(){
  cout << "--Begin " << __FUNCTION__ << endl;

  // Database for GEM
  TString db_template = fConfig->GetDBTemplate();
  
  TString input_name = fConfig->GetInputName();
  Ssiz_t first_t = input_name.Last('/') +1; // if a slash is not there, return 0.
  Ssiz_t last_t = input_name.Last('.');
  Int_t length_t = last_t - first_t;
  TString prefix_t = input_name(first_t,length_t);

  TString db_filename = Form("DBfiles/db_sbs.gems.dat_%s",prefix_t.Data());
  
  // Make a copy of template
  std::ifstream srce( db_template.Data(), std::ios::binary ) ;
  std::ofstream dest( db_filename.Data(), std::ios::binary ) ;
  dest << srce.rdbuf() ;
  srce.close();
  dest.close();
  cout << "--"<<__FUNCTION__ << ": DB template copied " << endl;
  
  cout << "--Opening  DB file: " << db_filename << endl;  
  FILE *db_file = fopen(db_filename.Data(),"a");
  
  // Insert a comment line 
  fprintf(db_file,"\n");
  fprintf(db_file,"# Pedestal Database for SLAC Beam Test %s",prefix_t.Data());
  
  // * Pre-Analysis GEM : Get Pedestal From Gaussian Fit

  const int nproj = 4;
  TH1D* hped_mean[nproj];
  TH1D* hped_rms[nproj];

  TString strProj[nproj]={"x1","y1","x2","y2"};
  Int_t sizeArray[nproj]={256, 512, 256, 512}; 

  // Summary histogram
  for(int iproj=0; iproj<nproj;iproj++){
      hped_mean[iproj] = new TH1D(Form("hped_mean_%d",iproj),
				  Form("Pedestal Mean vs strip, projection %s",strProj[iproj].Data()),
				  sizeArray[iproj],-0.5,sizeArray[iproj]-0.5);
      hped_rms[iproj] = new TH1D(Form("hped_rms_%d",iproj),
				 Form("Pedestal RMS vs strip, projection %s",strProj[iproj].Data()),
				 sizeArray[iproj],-0.5,sizeArray[iproj]-0.5);
  }

  Double_t ped_mean; //averaged by 6
  Double_t ped_rms; // averaged by sqrt(6)
  
  TString draw_text, hist_name;
  TH1D* h_fit = new TH1D("h_fit","Buffer historgram for pedestal fit",1e3,0.0,1e4);
  
  // Retrieve Channel Mapping from rootfiles ........
  double gem1_xstrip_id[256];
  double gem1_ystrip_id[512];
  double gem2_xstrip_id[256];
  double gem2_ystrip_id[512];
  double* strip_id[nproj] = { gem1_xstrip_id, gem1_ystrip_id,
			      gem2_xstrip_id, gem2_ystrip_id };

  // strip_id type should be an integer, but ....it was output as a float
  // Hardcoded by hand ,may need to think of a better solution to do this.
  // Anyway, I need pointers to load strip map
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  tree_raw->SetBranchAddress("sbs.gems.x1.strip",gem1_xstrip_id);
  tree_raw->SetBranchAddress("sbs.gems.y1.strip",gem1_ystrip_id);
  tree_raw->SetBranchAddress("sbs.gems.x2.strip",gem2_xstrip_id);
  tree_raw->SetBranchAddress("sbs.gems.y2.strip",gem2_ystrip_id);
  
  tree_raw->GetEntry(1); // in order to load strip map to these array
  
  // !!NOTE !! :
  // This routine will not work  if ZeroSuppression is TRUE in SBS-offline

  cout << "--Calculating Pedestals... " << endl;
  
  for(int iproj =0; iproj< nproj; iproj++){

    fprintf(db_file,"\nsbs.gems.%s.ped =",strProj[iproj].Data());
    int nch = sizeArray[iproj];

    for(int ich=0; ich<nch;ich++){

      // Get Strip id from array
      int strip = strip_id[iproj][ich]; 

      draw_text = Form("sbs.gems.%s.adc_sum[%d]",
		       strProj[iproj].Data(),
		       ich); // channel number
      tree_raw->Draw(Form("%s >> h_fit", draw_text.Data()),"","goff");
      
      GaussianFit(h_fit, ped_mean, ped_rms,iproj,strip);
	
      // File-Print pedestals
      if(ich%8==0)
	fprintf(db_file, " \\ \n");

      fprintf(db_file,"%d %.2f ",(int)strip, ped_mean);
      
      hped_mean[iproj]->SetBinContent(strip+1,ped_mean);
      hped_rms[iproj]->SetBinContent(strip+1,ped_rms);
      
      // counts++;
      // if(counts%10==0){
      // 	printf("\r Running %.1f %%  ",counts/total_counts*100);
      // }
      
    } 
  }
  
  cout <<"--Pedestal calibration is done !" <<endl;

  fclose(db_file);
  //Write  objects
  for(int iproj=0; iproj<nproj; iproj++){
    hped_mean[iproj]->Write();
    hped_rms[iproj]->Write();
  }

  return 0;
}

int BeamAnalysis::CalculateRMS(){
  // Header for GEM RMS

  TString input_name = fConfig->GetInputName();
  Ssiz_t first_t = input_name.Last('/') +1; // if a slash is not there, return 0.
  Ssiz_t last_t = input_name.Last('.');
  Int_t length_t = last_t - first_t;
  TString prefix_t = input_name(first_t,length_t);

  TString header_filename = Form("DBfiles/table_rms.h_%s",prefix_t.Data());
  FILE *header_file = fopen(header_filename.Data(),"w");

  fprintf(header_file,"\n");
  fprintf(header_file," //// RMS arrays run %s", prefix_t.Data());
  
  const int nproj = 4;
  TH1D* hped_mean[nproj];
  TH1D* hped_rms[nproj];

  TString strProj[nproj]={"x1","y1","x2","y2"};
  Int_t sizeArray[nproj]={256, 512, 256, 512}; 


  for(int iproj=0; iproj<nproj;iproj++){
    int nbins =sizeArray[iproj];
    hped_mean[iproj] = new TH1D(Form("hped_mean_%d",iproj),
				Form("Pedestal Mean vs strip,  %s",strProj[iproj].Data()),
				nbins,-0.5,nbins-0.5);
    hped_rms[iproj] = new TH1D(Form("hped_rms_%d",iproj),
			       Form("Pedestal RMS vs strip,  %s",strProj[iproj].Data()),
			       nbins,-0.5,nbins-0.5);
  }

  // Retrieve Channel Mapping from rootfiles ........
  double gem1_xstrip_id[256];
  double gem1_ystrip_id[512];
  double gem2_xstrip_id[256];
  double gem2_ystrip_id[512];
  double* strip_id[nproj] = { gem1_xstrip_id, gem1_ystrip_id,
			      gem2_xstrip_id, gem2_ystrip_id }; 
  // strip_id type should be an integer, but ....it was output as a float
  // Hardcoded by hand ,may need to think of a better solution to do this.
  // Anyway, I need a pointer to load strip map
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  tree_raw->SetBranchAddress("sbs.gems.x1.strip",gem1_xstrip_id);
  tree_raw->SetBranchAddress("sbs.gems.y1.strip",gem1_ystrip_id);
  tree_raw->SetBranchAddress("sbs.gems.x2.strip",gem2_xstrip_id);
  tree_raw->SetBranchAddress("sbs.gems.y2.strip",gem2_ystrip_id);
  
  tree_raw->GetEntry(1); // in order to load strip map to these array
  // Caution: And this needs to be fixed if ZeroSuppression is on in SBS-offline
  
  Double_t ped_mean; //averaged by 6
  Double_t ped_rms; // averaged by sqrt(6)
  
  TString draw_text, hist_name;
  TH1D* h_fit = new TH1D("h_fit","Buffer historgram for pedestal fit",5e2,-2e3,2e3);

  cout << "Calculating rms... " << endl;

  for(int iproj=0; iproj<4; iproj++){

    fprintf(header_file,"\n");
    fprintf(header_file,"double rms_%s[%d]={",strProj[iproj].Data(),sizeArray[iproj]);

    int nch = sizeArray[iproj];
    for(int ich=0; ich<nch;ich++){
      draw_text = Form("sbs.gems.%s.adc_sum[%d]-6*sbs.gems.%s.common_mode[%d]",
		       strProj[iproj].Data(),ich,
		       strProj[iproj].Data(),ich);
      tree_raw->Draw(Form("%s >> h_fit", draw_text.Data()),"","goff");
      h_fit->SetTitle(Form("Corrected Pedestal, %s, strip %d",
			   strProj[iproj].Data(), ich));
      
      GaussianFit(h_fit, ped_mean, ped_rms, iproj,ich);
      
      if(ich%16==0)
	fprintf(header_file,"\n");
      fprintf(header_file,"%.2f",ped_rms);
      if(ich!=nch-1)
	fprintf(header_file,", ");
      
      int strip = strip_id[iproj][ich];
      hped_mean[iproj]->SetBinContent(strip+1,ped_mean);
      hped_rms[iproj]->SetBinContent(strip+1,ped_rms);
    }
    
    fprintf(header_file,"};\n");
  }

  printf("RMS Calculation is Done ! \n");

  return 0;
}

int BeamAnalysis::Analysis(Bool_t kPlot){
  if (!kPlot)
    cout << "No Plot" <<endl;
  else
    cout << "Yes Plot" << endl;
  return 0;
}


void BeamAnalysis::GaussianFit(TH1D *h_fit, Double_t &mean, Double_t &sigma,
			       int iproj, int strip){
  
  Int_t bin_max = h_fit->GetMaximumBin();
  Double_t bincenter = h_fit->GetBinCenter(bin_max);
  Double_t bin_content_max = h_fit->GetBinContent(bin_max);
  Double_t rms = 100.0; // An initial guess

  Double_t par[3]; 
  par[0] = bin_content_max;
  par[1] = bincenter;
  par[2] = rms; 

  TF1 *f_gaus = new TF1("f_gaus","gaus",0,5e4);
  f_gaus->SetParameters(par);
  h_fit->Fit("f_gaus","QNR","",bincenter-2*rms,bincenter+2*rms);

  mean = f_gaus->GetParameter(1);
  sigma  = f_gaus->GetParameter(2);
  
  h_fit->Fit("f_gaus","QNR","",mean-2*sigma,mean+2*sigma);

  mean = f_gaus->GetParameter(1)/6.0; // averaged by 6
  sigma  = f_gaus->GetParameter(2)/sqrt(6);  // averaged by sqrt(6), assuming samples are not correlated
  // TCanvas *c1  = new TCanvas("c1","c1",800,800);
  // c1->cd();
  // h_fit->Draw();
  // f_gaus->Draw("same");
  // c1->SaveAs(Form("plots/FitPed-%d-%d.png",iproj,strip));
  // delete c1;
}

void BeamAnalysis::PrintSummary(){

}
