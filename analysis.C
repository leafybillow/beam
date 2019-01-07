
#define N_GEM 2; // Number of GEMS
#define THRESHOLD 3.0; // threshold for signal, e.g. 3 means 3*sigma 
#define N_CH_CMN 100;  // Number of strips used for common mode noise calculation
#define PITCH 400.0;   // pitch width, unit: um
#define CUT_MPL 1;     // cut on multiplicity, reject single wire hits (most likely noises or cross-talks)
#define CUT_COR 1.0;     // cut on correlation, based on correlation studies

Int_t analysis(TString filename="test.root", Bool_t isDebug = 0){
  
  gSystem->Load("libbeam.so");

  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "Raw ROOTFile " << filename << " is opened. " << endl;

  TString ped_filename = filename.ReplaceAll(".root","_ped.root");
  TFile* rf_ped = TFile::Open(ped_filename.Data());
  cout << "Pedestal ROOTFile " << ped_filename << " is opened. " << endl;

  TString output_filename = filename.ReplaceAll("_ped","_analyzed");
  TFile* rf_rec = TFile::Open(output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;

  // Create Trees and Branches
  TTree* tree_rec = new TTree("Rec","Rec"); // Reconstructed Tree
  
  BeamGEMPlane* gemPlane1 = new BeamGEMPlane();
  BeamGEMPlane* gemPlane2 = new BeamGEMPlane();
  BeamGEMTracker* tracker = new BeamGEMTracker();
  struct QDC{
    double us_lo;
    double us_hi;
    double ds_lo;
    double ds_hi;
  }qdc;

  TBranch* branch_qdc = tree_rec->Branch("qdc",&qdc,"us_lo/D:us_hi:ds_lo:ds_hi");
  TBranch* branch_gem1 = tree_rec->Branch("gem1",&gemPlane1);
  TBranch* branch_gem2 = tree_rec->Branch("gem2",&gemPlane2);
  TBranch* branch_track = tree_rec->Branch("tracker",&tracker);
  //_________________________________________________________________________________________

  // GEM common Variables
  const int ngem = N_GEM;
  const int napv = 6;  // number of apvs for each GEM
  const int nch = 128; // number of channel on each APV cards

  TString apv_tag[napv] = {"x","x",
			   "y","y","y","y"};
  Int_t base[napv] = {0, 128, 
		      0, 128 , 256, 384};
  Double_t cmn[ngem][napv];
  Int_t nch_cmn = N_CH_CMN; // number of channels used for computing  common mode

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

    //** Retrieve replayed data from the raw tree
    for(int igem=0; igem<2; igem++){
      for(int iapv=0; iapv<napv; iapv++){
	for(int ich=0; ich<nch; ich++){
	  
	}
      }
    }
    //** Calculate Common Mode 

    //** Do Pedestal and Common Mode correction 
    //** Calculate amplitude from each stripes and fill into a 1D histogram
    //** SearchHits in each readout coordinate from the 1D histogram
    //*** TestSplit 
    //** Verify  nhits match betweeen two readout coordinate
    //** Calculate cluster charge in 1D
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
  rf_ped->Close();
  tree_rec->Write();
  rf_rec->Close();
  return 0;
}

// User Define Functions

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


