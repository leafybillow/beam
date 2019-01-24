#include "analysis_rms.h"
#define THRESHOLD 3.0; // threshold for signal, e.g. 3 means 3*sigma 

Int_t analysis(TString filename="test.root", Bool_t isDebug = 0, Bool_t zeroSuppression=1){
  
  gSystem->Load("libbeam.so");

  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "Raw ROOTFile " << filename << " is opened. " << endl;

  Ssiz_t first_t = filename.Last('/') +1; // if a slash is not there, return 0.
  Ssiz_t last_t = filename.Last('.');
  Int_t length_t = last_t - first_t;
  TString prefix_t = filename(first_t,length_t);

  if(isDebug!=1){
    TString output_filename = prefix_t + "_analyzed.root";
    TFile* rf_rec = TFile::Open("rootfiles/"+output_filename,"RECREATE");
    cout << "ROOTFile " << output_filename << " is recreated. " << endl;
  }
  else{
    cout << "Debug mode is on " <<endl;
  }
  //__________________________________________________________________________________
  // New Reconstructed Tree
  TTree* tree_rec = new TTree("Rec","Rec"); 
  //And Build  Reconstruction Branch

  vector <double> vCharge_x1, vCharge_y1;
  vector <double> vPosition_x1, vPosition_y1;
  vector <double> vWidth_x1, vWidth_y1;
  vector <int> vSplit_x1, vSplit_y1;

  vector <double> vCharge_x2, vCharge_y2;
  vector <double> vPosition_x2, vPosition_y2;
  vector <double> vWidth_x2, vWidth_y2;
  vector <int> vSplit_x2, vSplit_y2;

  Int_t nHits_1, nHits_2;

  struct QDC{
    double us_lo; // upstream detector , low range, hi-sensitivity
    double us_hi; // upstream detector, high range,
    double ds_lo; // downstream detector
    double ds_hi; // downstream detector
  }qdc;

  TBranch* branch_qdc = tree_rec->Branch("qdc",&qdc,"us_lo/D:us_hi:ds_lo:ds_hi");

  tree_rec->Branch("gem1.charge_x",&vCharge_x1);
  tree_rec->Branch("gem1.charge_y",&vCharge_y1);
  tree_rec->Branch("gem1.position_x",&vPosition_x1);
  tree_rec->Branch("gem1.position_y",&vPosition_y1);
  tree_rec->Branch("gem1.width_x",&vWidth_x1);
  tree_rec->Branch("gem1.width_y",&vWidth_y1);
  tree_rec->Branch("gem1.split_x",&vSplit_x1);
  tree_rec->Branch("gem1.split_y",&vSplit_y1);
  tree_rec->Branch("gem1.nHits",&nHits_1);

  tree_rec->Branch("gem2.charge_x",&vCharge_x2);
  tree_rec->Branch("gem2.charge_y",&vCharge_y2);
  tree_rec->Branch("gem2.position_x",&vPosition_x2);
  tree_rec->Branch("gem2.position_y",&vPosition_y2);
  tree_rec->Branch("gem2.width_x",&vWidth_x2);
  tree_rec->Branch("gem2.width_y",&vWidth_y2);
  tree_rec->Branch("gem2.split_x",&vSplit_x2);
  tree_rec->Branch("gem2.split_y",&vSplit_y2);
  tree_rec->Branch("gem2.nHits",&nHits_2);

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
	// effectively zero suppression
	if(!zeroSuppression || bgStrip->GetADCsum()>THRESHOLD*sqrt(6)*rms[iProj][ich]){ 
	  bgStrip->Process();
	  bgProjection[iProj]->AddStrip(bgStrip);
	}
      } // End channel loop
      if(!isDebug){
	bgProjection[iProj]->Process();
      }
    } // End Projection Loop
    BeamGEMTracker* bgTracker = new BeamGEMTracker();
    BeamGEMPlane* bgPlane1 = new BeamGEMPlane("gem1");
    BeamGEMPlane* bgPlane2 = new BeamGEMPlane("gem2");
    bgPlane1->AddProjectionX(bgProjection[0]);
    bgPlane1->AddProjectionY(bgProjection[1]);
    bgPlane2->AddProjectionX(bgProjection[2]);
    bgPlane2->AddProjectionY(bgProjection[3]);

    if(!isDebug){
      bgPlane1->Process();
      bgPlane2->Process();
    }
    
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

    
    if (!isDebug){
      tree_rec->Fill();
    }

    if(isDebug){
	bgTracker->PlotResults(prefix_t,ievt);
    }

  } // End Event loop

  //__________________________________________________________________________________
  
  rf_raw->Close();
  if(isDebug!=1){
    tree_rec->Write();
    rf_rec->Close();
  }

  return 0;
}
