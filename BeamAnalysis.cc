#include "BeamGEMData.h"
#include "BeamGEMStrip.h"
#include "BeamGEMPlane.h"
#include "BeamGEMProjection.h"
#include "BeamGEMTracker.h"
#include "BeamAnalysis.h"
#include "BeamConfig.h"

#include "table_rms.h"

#include "TF1.h"
#include "TTree.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>

#define THRESHOLD 3.0

using namespace std;

ClassImp(BeamAnalysis);

BeamAnalysis::BeamAnalysis(BeamConfig *beamConfig){
  fConfig = beamConfig;

  kPlot=fConfig->GetPlotMode();
  
  rf_raw = TFile::Open(fConfig->GetInputName());
  TString input_filename = rf_raw->GetName(); 
  cout << "--Input file " << input_filename << " is opened. " << endl;
  if(!kPlot){
    rf_output = TFile::Open(fConfig->GetOutputName() ,"RECREATE");
    TString output_filename = rf_output->GetName();
    cout << "--Output ROOTFile " << output_filename << " is recreated. " << endl;
  }
  else
    cout << "--Plot Mode is ON" << endl;
}

BeamAnalysis::~BeamAnalysis(){
}

int BeamAnalysis::Process(){
  
  anaType = fConfig->GetAnalysisType();

  switch(anaType){
  case 0:
    Analysis();
    break;
  case 1:
    CalculatePed();
    break;
  case 2:
    CalculateRMS();
    break;
  }
    
  rf_raw->Close();
  if(!kPlot)
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

int BeamAnalysis::Analysis(){
  
  TString input_name = fConfig->GetInputName();
  Ssiz_t first_t = input_name.Last('/') +1; // if a slash is not there, return 0.
  Ssiz_t last_t = input_name.Last('.');
  Int_t length_t = last_t - first_t;
  TString prefix_t = input_name(first_t,length_t);

  cout << "--Begin Reconstruction and Tracking analysis " << endl;
  TTree *tree_rec;

  vector <double> vCharge_x1, vCharge_y1;
  vector <double> vPosition_x1, vPosition_y1;
  vector <double> vWidth_x1, vWidth_y1;
  vector <int> vSplit_x1, vSplit_y1;

  vector <double> vCharge_x2, vCharge_y2;
  vector <double> vPosition_x2, vPosition_y2;
  vector <double> vWidth_x2, vWidth_y2;
  vector <int> vSplit_x2, vSplit_y2;

  double gem1_baseline_rms[2];
  double gem1_baseline_mean[2];
  
  double gem2_baseline_rms[2];
  double gem2_baseline_mean[2];  // [0]: projX; [1]:projY
  
  Int_t nHits_1, nHits_2;

  struct QDC{
    double us_lo; // upstream detector , low range, hi-sensitivity
    double us_hi; // upstream detector, high range,
    double ds_lo; // downstream detector
    double ds_hi; // downstream detector
  }qdc;
  
  if(!kPlot){    
    // Reconstructed Tree
    tree_rec = new TTree("Rec","Rec"); 
    //And Build  Reconstruction Branch
    tree_rec->Branch("qdc",&qdc,"us_lo/D:us_hi:ds_lo:ds_hi");

    tree_rec->Branch("gem1.charge_x",&vCharge_x1);
    tree_rec->Branch("gem1.charge_y",&vCharge_y1);
    tree_rec->Branch("gem1.position_x",&vPosition_x1);
    tree_rec->Branch("gem1.position_y",&vPosition_y1);
    tree_rec->Branch("gem1.width_x",&vWidth_x1);
    tree_rec->Branch("gem1.width_y",&vWidth_y1);
    tree_rec->Branch("gem1.split_x",&vSplit_x1);
    tree_rec->Branch("gem1.split_y",&vSplit_y1);
    tree_rec->Branch("gem1.nHits",&nHits_1);
    tree_rec->Branch("gem1.ped_mean",gem1_baseline_mean,
		     "gem1.ped_mean[2]/D");
    tree_rec->Branch("gem1.ped_rms",gem1_baseline_rms,
		     "gem1.ped_rms[2]/D");
    
    tree_rec->Branch("gem2.charge_x",&vCharge_x2);
    tree_rec->Branch("gem2.charge_y",&vCharge_y2);
    tree_rec->Branch("gem2.position_x",&vPosition_x2);
    tree_rec->Branch("gem2.position_y",&vPosition_y2);
    tree_rec->Branch("gem2.width_x",&vWidth_x2);
    tree_rec->Branch("gem2.width_y",&vWidth_y2);
    tree_rec->Branch("gem2.split_x",&vSplit_x2);
    tree_rec->Branch("gem2.split_y",&vSplit_y2);
    tree_rec->Branch("gem2.nHits",&nHits_2);
    tree_rec->Branch("gem2.ped_mean",gem2_baseline_mean,
		     "gem2.ped_mean[2]/D");
    tree_rec->Branch("gem2.ped_rms",gem2_baseline_rms,
		     "gem2.ped_rms[2]/D");
  }
  //__________________________________________________________________________________
  // GEM Configuration Parameters

  Int_t nProj = 4; // number of projections, 2 GEM *(X+Y) = 4
  Int_t nadc = 6; // number of adc samples

  TString strADC[6]={"adc0","adc1","adc2","adc3","adc4","adc5"};
  TString strProj[4]={"x1","y1","x2","y2"};
  Int_t sizeArray[4]={256, 512, 256, 512}; 
  //__________________________________________________________________________________
  // Initialize EventReader for Raw Tree
  Double_t* rms[4]={rms_x1,rms_y1,rms_x2,rms_y2};

  Double_t us_hi,us_lo,ds_hi,ds_lo; // dummy variables
  // Note: qdc channels need to be changed depending on run configuration
  // Here just randomly pick up two channels for demonstration
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  tree_raw->SetBranchAddress("sbs.sbuscint.hadc2",&us_hi);
  tree_raw->SetBranchAddress("sbs.sbuscint.ladc2",&us_lo);
  tree_raw->SetBranchAddress("sbs.sbuscint.hadc1",&ds_hi);
  tree_raw->SetBranchAddress("sbs.sbuscint.ladc1",&ds_lo);

  BeamGEMData bgData[4];  // GEM Data Containers

  for(Int_t iProj=0;iProj<nProj;iProj++){

    bgData[iProj] = BeamGEMData(iProj,sizeArray[iProj]);

    for(Int_t iadc=0;iadc<nadc;iadc++){
      TString strBranch = Form("sbs.gems.%s.%s",
			       strProj[iProj].Data(),
			       strADC[iadc].Data());
      tree_raw->SetBranchAddress(strBranch, bgData[iProj].adc[iadc]);
    }
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.strip",strProj[iProj].Data()),
			       bgData[iProj].id_strip);
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.nch",strProj[iProj].Data()),
			       &(bgData[iProj].nChannel));
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.adc_sum",strProj[iProj].Data()),
			       bgData[iProj].adc_sum);
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.common_mode",strProj[iProj].Data()),
			       bgData[iProj].common_mode);
  }
  //__________________________________________________________________________________

  //*Event loop, reconstruction
  Bool_t kZeroSuppression = 1; // Suppressed by defaultp
  BeamGEMProjection* bgProjection[4];

  Int_t nentries = tree_raw->GetEntries();
  for(Int_t ievt=0;ievt<nentries;ievt++){
    if(ievt%200==0)
      cout << ievt << " events analyzed"  << endl;
    tree_raw->GetEntry(ievt);

    //** Retrieve replayed data from the raw tree
    //*** QDC
    qdc.us_lo = us_lo;
    qdc.us_hi = us_hi;
    qdc.ds_lo = ds_lo;
    qdc.ds_hi = ds_hi;

    //*** GEM
    for(Int_t iProj=0;iProj<nProj;iProj++){
      Int_t nChannel = (Int_t)bgData[iProj].nChannel;

      bgProjection[iProj] = new BeamGEMProjection( strProj[ bgData[iProj].id_Proj ],
						   (Int_t)bgData[iProj].nChannel);
      for(Int_t ich=0; ich<nChannel;ich++){
	  
	Double_t arADC[6];  	  // *** A Buffer container ADC Samples
	// common mode correction here
	for(Int_t iadc=0 ;iadc<nadc;iadc++){
	  arADC[iadc] = bgData[iProj].adc[iadc][ich] - bgData[iProj].common_mode[ich];
	} 
	Int_t myStripID = (Int_t)bgData[iProj].id_strip[ich];
	BeamGEMStrip* bgStrip = new BeamGEMStrip(arADC,myStripID);
	// Zero suppression
	
	if(bgStrip->GetADCsum()>THRESHOLD*sqrt(6)*rms[iProj][ich])
	  kZeroSuppression = 0; // Not Suppressed
	else
	  kZeroSuppression = 1;
	
	bgStrip->SetZeroSuppression(kZeroSuppression);

	// bgStrip->Process(); // FIXME: Not implemented now
	bgProjection[iProj]->AddStrip(bgStrip);

      } // End channel loop
      bgProjection[iProj]->Process();
    } // End Projection Loop
    BeamGEMTracker* bgTracker = new BeamGEMTracker();
    BeamGEMPlane* bgPlane1 = new BeamGEMPlane("gem1");
    BeamGEMPlane* bgPlane2 = new BeamGEMPlane("gem2");
    bgPlane1->AddProjectionX(bgProjection[0]);
    bgPlane1->AddProjectionY(bgProjection[1]);
    bgPlane2->AddProjectionX(bgProjection[2]);
    bgPlane2->AddProjectionY(bgProjection[3]);

    bgPlane1->Process();
    bgPlane2->Process();
    
    bgTracker->AddPlane(bgPlane1);
    bgTracker->AddPlane(bgPlane2);
    // FIXME: a better way to do this

    vCharge_x1 = bgPlane1->GetChargeX();
    vCharge_y1 = bgPlane1->GetChargeY();
    vPosition_x1 = bgPlane1->GetPositionX();
    vPosition_y1 = bgPlane1->GetPositionY();

    vWidth_x1 = bgPlane1->GetWidthX();
    vWidth_y1 = bgPlane1->GetWidthY();

    nHits_1 = bgPlane1->GetNHits();

    vCharge_x2 = bgPlane2->GetChargeX();
    vCharge_y2 = bgPlane2->GetChargeY();
    vPosition_x2 = bgPlane2->GetPositionX();
    vPosition_y2 = bgPlane2->GetPositionY();

    vWidth_x2 = bgPlane2->GetWidthX();
    vWidth_y2 = bgPlane2->GetWidthY();

    nHits_2 = bgPlane2->GetNHits();

    gem1_baseline_mean[0]
      = bgPlane1->GetProjectionX()->GetBaselineMean();
    gem1_baseline_mean[1]
      = bgPlane1->GetProjectionY()->GetBaselineMean();

    gem1_baseline_rms[0]
      = bgPlane1->GetProjectionX()->GetBaselineRMS();
    gem1_baseline_rms[1]
      = bgPlane1->GetProjectionY()->GetBaselineRMS();

    gem2_baseline_mean[0]
      = bgPlane2->GetProjectionX()->GetBaselineMean();
    gem2_baseline_mean[1]
      = bgPlane2->GetProjectionY()->GetBaselineMean();

    gem2_baseline_rms[0]
      = bgPlane2->GetProjectionX()->GetBaselineRMS();
    gem2_baseline_rms[1]
      = bgPlane2->GetProjectionY()->GetBaselineRMS();

    if(!kPlot){
      tree_rec->Fill();
    }

    if(kPlot){
      if(gem1_baseline_rms[1]>100 )
      // bgPlane1->GetProjectionY()->PlotResults(prefix_t,ievt);
      bgTracker->PlotResults(prefix_t,ievt);
    }

  } // End Event loop
  
  if(!kPlot)
    tree_rec->Write();

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
  if(kPlot){
    TCanvas *c1  = new TCanvas("c1","c1",800,800);
    c1->cd();
    h_fit->Draw();
    f_gaus->Draw("same");
    c1->SaveAs(Form("plots/FitPed-%d-%d.png",iproj,strip));
    delete c1;
  }
}

void BeamAnalysis::PrintSummary(){

}
