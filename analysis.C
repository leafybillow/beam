#include "analysis.h"

#define N_GEM 2; // Number of GEMS
#define THRESHOLD 3.0; // threshold for signal, e.g. 3 means 3*sigma 
#define N_CH_CMN 100;  // Number of strips used for common mode noise calculation
#define PITCH 400.0;   // pitch width, unit: um
#define CUT_MPL 1;     // cut on multiplicity, reject single wire hits (most likely noises or cross-talks)
#define CUT_COR 1.0;     // cut on correlation, based on correlation studies

Int_t analysis(TString filename="test.root", Bool_t isDebug = 0){
  
  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "ROOTFile " << filename << " is opened. " << endl;

  TString output_filename = filename.ReplaceAll(".root","_analyzed.root");
  TFile* rf_rec = TFile::Open(output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;

  // Create Trees and Branches
  TTree *tree_rec = new TTree("tree_rec","Ped"); // Reconstructed Tree
  
  GEMPlane gemPlane_1, gemPlane_2;
  Track track;
  QDC qdc;
  leaflist_track.ReplaceAll("N_HITS_MAX",Form("%d",N_HITS_MAX));
  leaflist_gem.ReplaceAll("N_HITS_MAX",Form("%d",N_HITS_MAX));
  // cout << leaflist_gem << endl;
  // cout << leaflist_track << endl;
  TBranch *branch_qdc = tree_rec->Branch("qdc",&qdc,leaflist_qdc);
  TBranch *branch_gem1 = tree_rec->Branch("gem1",&gemPlane_1,leaflist_gem);
  TBranch *branch_gem2 = tree_rec->Branch("gem2",&gemPlane_2,leaflist_gem);
  TBranch *branch_track = tree_rec->Branch("track",&track,leaflist_track);
  //_________________________________________________________________________________________

  // * Pre-Analysis GEM : Get Pedestal From Gaussian Fit
  const int ngem = N_GEM;
  const int napv = 6;  // number of apvs for each GEM
  const int nch = 128; // number of channel on each APV cards

  TString apv_tag[napv] = {"x","x",
			   "y","y","y","y"};
  Int_t base[napv] = {0, 128, 
		      0, 128 , 256, 384};
  Double_t ped_mean[ngem][napv][nch]; // [napv][nch] : averaged by 6
  Double_t ped_rms[ngem][napv][nch]; //[napv][nch]: averaged by sqrt(6)
  
  TString draw_text, hist_name;
  TH1D* h_ped = new TH1D("h_ped","Buffer historgram for pedestal fit",1e4,-0.5,2e4-0.5);

  for(int igem =0; igem< ngem; igem++){
    for(int iapv=0; iapv< napv; iapv++){
      for(int ich=0; ich<nch;ich++){
	draw_text = Form("sbs.gems.%s%d.adc_sum[%d]",
			 apv_tag[iapv].Data(),
			 igem+1,
			 ich+base[iapv]); // stripe number
	// tree_raw->Draw(Form("%s >> h_ped", draw_text.Data()),"","goff");
	// PedestalFit(h_ped,
	// 	    ped_mean[igem][iapv][ich],
	// 	    ped_rms[igem][iapv][ich]);
	// cout << igem << "\t" << iapv << "\t" << ich << "\t" ;
	// cout << ped_mean[igem][iapv][ich]  << "\t" << ped_rms[igem][iapv][ich]<< endl;
      } 
    } 
  }
  // * End Pre-Analysis
  //_________________________________________________________________________________________
  
  Double_t us_hi,us_lo,ds_hi,ds_lo; // dummy variables
  tree_raw->SetBranchAddress("sbs.sbuscint.hadc0",&us_hi);
  tree_raw->SetBranchAddress("sbs.sbuscint.ladc0",&us_lo);
  tree_raw->SetBranchAddress("sbs.sbuscint.hadc1",&ds_hi);
  tree_raw->SetBranchAddress("sbs.sbuscint.ladc1",&ds_lo);
  //*Event loop, reconstruction
  Int_t nentries = tree_raw->GetEntries();
  for(int ievt=0;ievt<nentries;ievt++){
    tree_raw->GetEntry(ievt);
    // -- Analyze Detector/QDC signal
    //** Get QDC from both low and high range channels, and upstream and downstream detectors
    qdc.us_lo = us_lo;
    qdc.us_hi = us_hi;
    qdc.ds_lo = ds_lo;
    qdc.ds_hi = ds_hi;
    // -- Now, Analyze GEM 
    //** Calculate Common Mode 

    //** Retrieve replayed data from the raw tree
    //** Do Pedestal and Common Mode correction and fill 1D histograms
    //** SearchHits in each readout coordinate from histogram
    //*** TestSplit 
    //** Verify  nhits match betweeen two readout coordinate
    //** Calculate deposited charge in 1D
    //** Calculate Multiplicity in 1D 
    //** Calculate Position in 1D , using centroid for now
    //** Reconstruct 2D hits from correlation
    //** Build Tracks
    //*** Sort second GEM hits array to match to first GEM.
    //*** Calculate slopes
    //*** Calculate angles
    //*** Option: Plot Tracks

    //** Fill this entry
    tree_rec->Fill();
    //** Option: output plots
  }  // End Event loop
  //_________________________________________________________________________________________
  rf_raw->Close();
  
  tree_rec->Write();
  rf_rec->Close();
  return 0;
}

// User Define Functions
void PedestalFit(TH1D *h_ped, Double_t &mean, Double_t &sigma){

  Int_t bin_max = h_ped->GetMaximumBin();
  Double_t bincenter = h_ped->GetBinCenter(bin_max);
  Double_t bin_content_max = h_ped->GetBinContent(bin_max);
  Double_t rms = 100.0; // An initial guess

  Double_t par[3]; 
  par[0] = bin_content_max;
  par[1] = bincenter;
  par[2] = rms; 

  TF1 *f_gaus = new TF1("f_gaus","gaus",0,5e4);
  f_gaus->SetParameters(par);
  h_ped->Fit("f_gaus","QN","",bincenter-2*rms,bincenter+2*rms);
  
  mean = f_gaus->GetParameter(1)/6.0; // averaged by 6
  sigma  = f_gaus->GetParameter(2)/sqrt(6);  // averaged by sqrt(6)

  // Crude Process
  // mean = bincenter;
  // sigma = rms;
}


Double_t GetCommonMode(Int_t *arr_integral){
}

Double_t GetStripeAmpl(Int_t *arr_sample, Double_t ped ){
  return ampl;
}

TH1D* GetHistX(){
  return hist_x;
}

TH1D* GetHistY(){
  return hist_y;
}


