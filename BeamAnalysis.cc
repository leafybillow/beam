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
#include "TString.h"

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
    rf_output = TFile::Open(fConfig->GetOutputName(),
			    "RECREATE");
    TString output_filename = rf_output->GetName();
    
    cout << "--Output ROOTFile "
	 << output_filename
	 << " is recreated. " << endl;
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
  TString db_path = fConfig->GetDBPath();
  TString db_filename = Form("%sdb_sbs.gems.dat_%s",
			     db_path.Data(),
			     prefix_t.Data());
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

    strProj[iproj] = projKey[iproj];
    
    if(strProj[iproj].Contains("x"))
      sizeArray[iproj] = 256;
    else if(strProj[iproj].Contains("y"))
      sizeArray[iproj] = 512;
    
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
  TString db_path = fConfig->GetDBPath();
  TString header_filename = Form("%s/%s_rms.table",
				 db_path.Data(),
				 prefix_t.Data());
  
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

  LoadRMS();
  
  cout << "--Begin Reconstruction and Tracking analysis " << endl;
  TTree *tree_rec;
  const Int_t nproj = projKey.size();
  vector <vector<double> > vCharge; //[iproj][ihits]
  vector <vector<double> > vPosition;
  vector <vector<double> > vWidth;
  vector <vector<int> > vSplit;
  vector <double> charge_sum;
  
  vector <int> vNhits;
  vector <int> vNhits_gem;
  vector <double> baseline_rms;
  vector <double> baseline_mean;

  // Reconstruction Detector Hits
  vector< vector<double> > vDet_x;
  vector< vector<double> > vDet_y;
  vector< vector<double> > vDet_theta;
  vector< vector<double> > vDet_phi;

  vector <double> dummy_vec_double;
  vector <int> dummy_vec_int;
  double dummy_double;
  int dummy_int;

  for(int iproj=0;iproj<nproj;iproj++){
    vCharge.push_back(dummy_vec_double);
    vPosition.push_back(dummy_vec_double);
    vWidth.push_back(dummy_vec_double);
    vSplit.push_back(dummy_vec_int);
    baseline_rms.push_back(dummy_double);
    baseline_mean.push_back(dummy_double);
    charge_sum.push_back(dummy_double);
    vNhits.push_back(dummy_int);
    if(iproj%2==0)
      vNhits_gem.push_back(dummy_int);
  }

  int ndets = fConfig->GetNDets();
  for(int idet=0;idet<ndets ;idet++){
    vDet_x.push_back(dummy_vec_double);
    vDet_y.push_back(dummy_vec_double);
    vDet_theta.push_back(dummy_vec_double);
    vDet_phi.push_back(dummy_vec_double);

  }
  int nTracks;
  bool isGoldenTrack;
  // Initialize EventReader for Raw Tree
  vector<Int_t> vec_qdc_ch = fConfig->GetQDCChannel();
  Int_t n_qdc_ch = vec_qdc_ch.size();
  Double_t *det_qdc_lr = new Double_t[n_qdc_ch];
  Double_t *det_qdc_hr = new Double_t[n_qdc_ch];

  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  // Reconstructed Tree
  tree_rec = new TTree("Rec","Rec");

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

  
  if(!kPlot){
    for(int idet=0; idet<ndets;idet++){
      tree_rec->Branch(Form("det%d_x",idet+1),&vDet_x[idet]);
      tree_rec->Branch(Form("det%d_y",idet+1),&vDet_y[idet]);
      tree_rec->Branch(Form("det%d_theta",idet+1),&vDet_theta[idet]);
      tree_rec->Branch(Form("det%d_phi",idet+1),&vDet_phi[idet]);
    }
    tree_rec->Branch("ntracks",&nTracks);
    tree_rec->Branch("isGoldenTrack",&isGoldenTrack);
  }
  
  if(!kPlot){    
    //And Build  Reconstruction Branch
    for(int iproj=0;iproj<nproj;iproj++){
      TString str_key = projKey[iproj];
      const char *key = str_key.Data();
      tree_rec->Branch(Form("charge_%s",key),&vCharge[iproj]);
      tree_rec->Branch(Form("position_%s",key),&vPosition[iproj]);
      tree_rec->Branch(Form("width_%s",key),&vWidth[iproj]);
      tree_rec->Branch(Form("split_%s",key),&vSplit[iproj]);
      tree_rec->Branch(Form("ped_mean_%s",key),&baseline_mean[iproj]);
      tree_rec->Branch(Form("ped_rms_%s",key),&baseline_rms[iproj]);
      tree_rec->Branch(Form("charge_sum_%s",key),&charge_sum[iproj]);
      tree_rec->Branch(Form("nHits_%s",key),&vNhits[iproj]);
    }

    for(int igem=0;igem<n_gem;igem++){
      tree_rec->Branch(Form("nHits_gem%d",igem+1),&vNhits_gem[igem]);
    }
      
  }
  //_________________________________________________________________________________
  // GEM Configuration Parameters
  Int_t nadc = 6; // number of adc samples
  TString strADC[6]={"adc0","adc1","adc2","adc3","adc4","adc5"};
  //_________________________________________________________________________________
  BeamGEMData bgData[nproj];  // GEM Data Containers
  Int_t sizeArray;
  
  for(Int_t iproj=0;iproj<nproj;iproj++){
    
    if(projKey[iproj].Contains("x"))
      sizeArray=256;
    else if(projKey[iproj].Contains("y"))
      sizeArray=512;

    const char* key = projKey[iproj].Data();
    bgData[iproj] = BeamGEMData(iproj,sizeArray);

    for(Int_t iadc=0;iadc<nadc;iadc++){
      TString strBranch = Form("sbs.gems.%s.%s",
			       key,strADC[iadc].Data());
      tree_raw->SetBranchAddress(strBranch, bgData[iproj].adc[iadc]);
    }
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.strip",key),
			       bgData[iproj].id_strip);
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.nch",key),
			       &(bgData[iproj].nChannel));
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.adc_sum",key),
			       bgData[iproj].adc_sum);
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.common_mode",key),
			       bgData[iproj].common_mode);
  }
  //__________________________________________________________________________________

  //*Event loop, reconstruction
  Bool_t kZeroSuppression = 1; // Suppressed by default
  
  Int_t nentries = tree_raw->GetEntries();
  for(Int_t ievt=0;ievt<nentries;ievt++){
    if(ievt%200==0)
      cout << ievt << " events analyzed"  << endl;
    //** Retrieve replayed data from the raw tree
    tree_raw->GetEntry(ievt);

    //*** GEM
    BeamGEMTracker* bgTracker = new BeamGEMTracker();
    bgTracker->LoadDetectorGeometry( fConfig );
    bgTracker->LoadQDC(det_qdc_hr[0]);
    BeamGEMPlane* bgPlane[n_gem];
    BeamGEMProjection* bgProjection[nproj];

    for(int igem=0;igem<n_gem;igem++){
      bgPlane[igem] = new BeamGEMPlane();
    }
    for(Int_t iproj=0;iproj<nproj;iproj++){
      Int_t nChannel = (Int_t)bgData[iproj].nChannel;

      bgProjection[iproj] = new BeamGEMProjection( projKey[iproj],
						   (Int_t)bgData[iproj].nChannel);
      for(Int_t ich=0; ich<nChannel;ich++){
	  
	Double_t arADC[6];  	  // *** A Buffer container ADC Samples
	// common mode correction here
	for(Int_t iadc=0 ;iadc<nadc;iadc++){
	  arADC[iadc] = bgData[iproj].adc[iadc][ich] - bgData[iproj].common_mode[ich];
	} 
	Int_t myStripID = (Int_t)bgData[iproj].id_strip[ich];
	BeamGEMStrip* bgStrip = new BeamGEMStrip(arADC,myStripID);
	// Zero suppression
	
	if(bgStrip->GetADCsum()>zs_threshold*sqrt(6)*rms[iproj][ich])
	  kZeroSuppression = 0; // Not Suppressed
	else
	  kZeroSuppression = 1;
	
	bgStrip->SetZeroSuppression(kZeroSuppression);

	// bgStrip->Process(); // FIXME: Not implemented now
	bgProjection[iproj]->AddStrip(bgStrip);

      } // End channel loop
      bgProjection[iproj]->Process();

      TString this_key = projKey[iproj];
      TString proj_type = this_key[0];
      Int_t gem_id = (this_key.Remove(0,1)).Atoi(); 
      if(proj_type=="x"){
	bgPlane[gem_id-1]->SetID(gem_id);
	bgPlane[gem_id-1]->AddProjectionX(bgProjection[iproj]);
      }
      else if(proj_type=="y")
	bgPlane[gem_id-1]->AddProjectionY(bgProjection[iproj]);
      
    } // End Projection Loop
    vector<Double_t> gem_z ;
    gem_z = fConfig->GetZ_GEM();
    for(int igem=0;igem<n_gem;igem++){
      bgPlane[igem]->SetPositionZ(gem_z[igem]);
      bgPlane[igem]->Process();
      bgTracker->AddPlane(bgPlane[igem]);
    }
    
    bgTracker->Process();
    
    // Loading analysis results to new branch pointer
    vDet_x = bgTracker->GetDetX();
    vDet_y = bgTracker->GetDetY();
    vDet_theta = bgTracker->GetDetTheta();
    vDet_phi = bgTracker->GetDetPhi();
    nTracks = bgTracker->GetNTracks();
    isGoldenTrack = bgTracker->IsGoldenTrack();
    
    for(int iproj=0;iproj<nproj;iproj++){
      TString this_key = projKey[iproj];
      TString proj_type = this_key[0];
      Int_t gem_id = (this_key.Remove(0,1)).Atoi(); 

      BeamGEMPlane *this_plane = bgPlane[gem_id-1];
      if(proj_type=="x"){
	vCharge[iproj] = this_plane->GetChargeX();
	vPosition[iproj] = this_plane->GetPositionX();
	vWidth[iproj] = this_plane->GetWidthX();
	baseline_mean[iproj] = this_plane->GetProjectionX()->GetBaselineMean();
	baseline_rms[iproj] = this_plane->GetProjectionX()->GetBaselineRMS();
	charge_sum[iproj] = this_plane->GetProjectionX()->GetChargeSum();
	vNhits[iproj] = this_plane->GetProjectionX()->GetNHits();
      }
      else if(proj_type=="y"){
	vCharge[iproj] = this_plane->GetChargeY();
	vPosition[iproj] = this_plane->GetPositionY();
	vWidth[iproj] = this_plane->GetWidthY();
	baseline_mean[iproj] = this_plane->GetProjectionY()->GetBaselineMean();
	baseline_rms[iproj] = this_plane->GetProjectionY()->GetBaselineRMS();
	charge_sum[iproj] = this_plane->GetProjectionY()->GetChargeSum();
	vNhits[iproj] = this_plane->GetProjectionY()->GetNHits();
      }
	  
    }
    for(int igem=0;igem<n_gem;igem++){
      vNhits_gem[igem] =bgPlane[igem]->GetNHits();
    }

    if(!kPlot){
      tree_rec->Fill();
    }

    if(kPlot){
      bgTracker->PlotResults(prefix_t,ievt);
    }

  } // End Event loop
  
  if(!kPlot)
    rf_output->WriteTObject(tree_rec);

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

  TF1 *f_gaus = new TF1("f_gaus","gaus",-5e4,5e4);
  f_gaus->SetParameters(par);
  h_fit->Fit("f_gaus","QNR","",bincenter-rms,bincenter+rms);

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
  TString rf_path = fConfig->GetRFPath();
  TString filename = Form("%s%s_rms.root",
			  rf_path.Data(),
			  prefix_t.Data());
  cout << "--Loading RMS" << endl;
  
  TFile *rms_rootfile = TFile::Open(filename);
  // If rms does not exist
  if(rms_rootfile->IsOpen())
    cout << "--Opened RMS rootfile" << endl;
  
  TH1D *hrms_buff;
  vector< Double_t> vector_rms;

  Int_t nproj = projKey.size();
  for(int iproj=0;iproj<nproj;iproj++){
    hrms_buff = (TH1D*)rms_rootfile->Get(Form("hped_rms_%d",iproj));

    Int_t nbins = hrms_buff->GetEntries();
    for(int ibin=0;ibin<nbins;ibin++){
      Double_t bin_content = hrms_buff->GetBinContent(ibin+1);
      vector_rms.push_back(bin_content);
    }
    rms.push_back(vector_rms);
    vector_rms.clear();
  }
  rms_rootfile->Close();
  cout << "--Number of Projections found: " << rms.size() << endl;
  cout << "--Closed RMS rootfile" << endl;
  return 0;
}

void BeamAnalysis::PrintSummary(){

}








