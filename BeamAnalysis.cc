#include "BeamGEMData.h"
#include "BeamGEMStrip.h"
#include "BeamGEMPlane.h"
#include "BeamGEMProjection.h"
#include "BeamGEMTracker.h"
#include "BeamAnalysis.h"
#include "BeamConfig.h"
#include "BeamParameters.h"

#include "TF1.h"
#include "TTree.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>

using namespace std;

ClassImp(BeamAnalysis);

BeamAnalysis::BeamAnalysis(BeamConfig *beamConfig){
  fConfig = beamConfig;

  kPlot=fConfig->GetPlotMode();
  n_gem = fConfig->GetNGEMs();

  for(int igem = 0;igem<n_gem;igem++){
    projKey.push_back(Form("x%d",igem+1));
    projKey.push_back(Form("y%d",igem+1));
  }
  
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
  const int nproj = projKey.size();
  // Summary histogram
  TH1D* hped_mean[nproj];
  TH1D* hped_rms[nproj];

  TString strProj[nproj];
  Int_t sizeArray[nproj];

  for(int iproj=0; iproj<nproj;iproj++){
    sizeArray[iproj] = ( (iproj%2)==0 ? 256 : 512);

    strProj[iproj] = projKey[iproj];
    
    hped_mean[iproj] = new TH1D(Form("hped_mean_%s",strProj[iproj].Data() ),
				Form("Pedestal Mean vs strip, projection %s",strProj[iproj].Data()),
				sizeArray[iproj],-0.5,sizeArray[iproj]-0.5);
    hped_rms[iproj] = new TH1D(Form("hped_rms_%s",strProj[iproj].Data() ),
			       Form("Pedestal RMS vs strip, projection %s",strProj[iproj].Data()),
			       sizeArray[iproj],-0.5,sizeArray[iproj]-0.5);
  }

  Double_t ped_mean; //averaged by 6
  Double_t ped_rms; // averaged by sqrt(6)
  
  TString draw_text, hist_name;
  TH1D* h_fit = new TH1D("h_fit","Buffer historgram for pedestal fit",1e3,0.0,1e4);
  
  // Retrieve Channel Mapping from rootfiles ........
  // strip_id type should be an integer, but ....it was output as a float
  // Hardcoded by hand ,may need to think of a better solution to do this.
  // Anyway, I need pointers to load strip map
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  double* strip_id[nproj];
  for(int iproj=0;iproj<nproj;iproj++){
    strip_id[iproj] = new double[ sizeArray[iproj] ];
    TString branch_name = Form("sbs.gems.%s.strip",strProj[iproj].Data());
    tree_raw->SetBranchAddress(branch_name,strip_id[iproj]);
  }
  
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

  TString header_filename = Form("DBfiles/run%s_rms.table",prefix_t.Data());
  FILE *header_file = fopen(header_filename.Data(),"w");

  fprintf(header_file,"\n");
  fprintf(header_file," ### RMS arrays run %s", prefix_t.Data());
  
  const int nproj = 2*n_gem;
  TH1D* hped_mean[nproj];
  TH1D* hped_rms[nproj];

  TString strProj[nproj];
  Int_t sizeArray[nproj];


  for(int iproj=0; iproj<nproj;iproj++){
    sizeArray[iproj] = ( (iproj%2)==0 ? 256 : 512);
    TString strDim = ( (iproj%2)==0 ? "x" : "y");
    strProj[iproj] = Form("%s%d",strDim.Data(), (iproj/2 +1));

    int nbins =sizeArray[iproj];
    hped_mean[iproj] = new TH1D(Form("hped_mean_%d",iproj),
				Form("Pedestal Mean vs strip,  %s",strProj[iproj].Data()),
				nbins,-0.5,nbins-0.5);
    hped_rms[iproj] = new TH1D(Form("hped_rms_%d",iproj),
			       Form("Pedestal RMS vs strip,  %s",strProj[iproj].Data()),
			       nbins,-0.5,nbins-0.5);
  }

  // Retrieve Channel Mapping from rootfiles ........
  // strip_id type should be an integer, but ....it was output as a float
  // Hardcoded by hand ,may need to think of a better solution to do this.
  // Anyway, I need a pointer to load strip map
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  double* strip_id[nproj];
  for(int iproj=0;iproj<nproj;iproj++){
    strip_id[iproj] = new double[ sizeArray[iproj] ];
    TString branch_name = Form("sbs.gems.%s.strip",strProj[iproj].Data());
    tree_raw->SetBranchAddress(branch_name,strip_id[iproj]);
  }
  
  tree_raw->GetEntry(1); // in order to load strip map to these array
  // Caution: And this needs to be fixed if ZeroSuppression is on in SBS-offline
  
  Double_t ped_mean; //averaged by 6
  Double_t ped_rms; // averaged by sqrt(6)
  
  TString draw_text, hist_name;
  TH1D* h_fit = new TH1D("h_fit","Buffer historgram for pedestal fit",5e2,-2e3,2e3);

  cout << "Calculating rms... " << endl;

  for(int iproj=0; iproj<nproj; iproj++){

    fprintf(header_file,"\n");
    fprintf(header_file,"rms_%s[%d]={",strProj[iproj].Data(),sizeArray[iproj]);

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
  
  //Write  objects
  for(int iproj=0; iproj<nproj; iproj++){
    hped_mean[iproj]->Write();
    hped_rms[iproj]->Write();
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
  Int_t nProj = 2*n_gem;
  vector <vector<double> > vCharge_x, vCharge_y;
  vector <vector<double> > vPosition_x, vPosition_y;
  vector <vector<double> > vWidth_x, vWidth_y;
  vector <vector<int> > vSplit_x, vSplit_y;
  
  vector <int> vNhits;

  vector <double> baseline_rms_x;
  vector <double> baseline_rms_y;
  vector <double> baseline_mean_x;
  vector <double> baseline_mean_y;
  
  vector <double> dummy_vec_double;
  vector <int> dummy_vec_int;
  double dummy_double;
  int dummy_int;

  for(int igem=0;igem<ngem;igem++){
    vCharge_x.push_back(dummy_vec_double);
    vCharge_y.push_back(dummy_vec_double);
    vPosition_x.push_back(dummy_vec_double);
    vPosition_y.push_back(dummy_vec_double);
    vWidth_x.push_back(dummy_vec_double);
    vWidth_y.push_back(dummy_vec_double);
    vSplit_x.push_back(dummy_vec_int);
    vSplit_y.push_back(dummy_vec_int);

    baseline_rms_x.push_back(dummy_double);
    baseline_rms_y.push_back(dummy_double);
    baseline_mean_x.push_back(dummy_double);
    baseline_mean_y.push_back(dummy_double);
    vNhits.push_back(dummy_int);
  }
  
  // Initialize EventReader for Raw Tree
  // Double_t* rms[4]={rms_x1,rms_y1,rms_x2,rms_y2};
  LoadRMS();

  vector<Int_t> vec_qdc_ch = fConfig->GetQDCChannel();
  Int_t n_qdc_ch = vec_qdc_ch.size();
  Double_t *det_qdc_lr = new Double_t[n_qdc_ch];
  Double_t *det_qdc_hr = new Double_t[n_qdc_ch];
  
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  // Reconstructed Tree
  tree_rec = new TTree("Rec","Rec");

  if(!kPlot){
    for(int iqdc=0;iqdc<n_qdc_ch;iqdc++){

      Int_t qdc_ch = vec_qdc_ch[iqdc];
      TString lrqdc_name = Form("sbs.sbuscint.ladc%d",qdc_ch);
      TString hrqdc_name = Form("sbs.sbuscint.hadc%d",qdc_ch);
      tree_raw->SetBranchAddress(lrqdc_name.Data(), &det_qdc_lr[iqdc]);
      tree_raw->SetBranchAddress(hrqdc_name.Data(), &det_qdc_hr[iqdc]);
      //   ||     ||
      // This is a bridge
      //   ||     ||
      tree_rec->Branch(Form("det%d_qdc_lr",iqdc+1),&det_qdc_lr[iqdc]);
      tree_rec->Branch(Form("det%d_qdc_hr",iqdc+1),&det_qdc_hr[iqdc]);
    }
  }  

  if(!kPlot){    
    //And Build  Reconstruction Branch
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


  Int_t nadc = 6; // number of adc samples

  TString strADC[6]={"adc0","adc1","adc2","adc3","adc4","adc5"};
  TString strProj[4]={"x1","y1","x2","y2"};
  Int_t sizeArray[4]={256, 512, 256, 512}; 
  //__________________________________________________________________________________
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
  Bool_t kZeroSuppression = 1; // Suppressed by default
  BeamGEMProjection* bgProjection[4];

  Int_t nentries = tree_raw->GetEntries();
  for(Int_t ievt=0;ievt<nentries;ievt++){
    if(ievt%200==0)
      cout << ievt << " events analyzed"  << endl;
    //** Retrieve replayed data from the raw tree
    tree_raw->GetEntry(ievt);

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
	
	if(bgStrip->GetADCsum()>zs_threshold*sqrt(6)*rms[iProj][ich])
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
  Double_t rms = 50.0; // An initial guess

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

int BeamAnalysis::LoadRMS(){
  TString input_name = fConfig->GetInputName();
  Ssiz_t first_t = input_name.Last('/') +1; // if a slash is not there, return 0.
  Ssiz_t last_t = input_name.Last('.');
  Int_t length_t = last_t - first_t;
  TString prefix_t = input_name(first_t,length_t);

  TString filename = Form("rootfiles/%s_rms.root",prefix_t.Data());

  TFile *rms_rootfile = TFile::Open(filename);
  // If rms does not exist
  
  TH1D *hrms_buff;
  vector< Double_t> vector_rms;

  Int_t nproj = projKey.size();
  for(int iproj=0;iproj<nproj;iproj++){
    hrms_buff = (TH1D*)rms_rootfile->FindObject(Form("hped_rms_%s",projKey[iproj].Data()));

    Int_t nbins = hrms_buff->GetNbinsX();
    for(int ibin=0;ibin<nbins;ibin++){
      Double_t bin_content = hrms_buff->GetBinContent(ibin+1);
      vector_rms.push_back(bin_content);
    }
    rms.push_back(vector_rms);
    vector_rms.clear();
  }
  rms_rootfile->Close();
  return 0;
}

void BeamAnalysis::PrintSummary(){

}
