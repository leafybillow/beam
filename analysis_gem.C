#include "analysis_rms.h"
#define THRESHOLD 3.0; // threshold for signal, e.g. 3 means 3*sigma 

Int_t analysis_gem(TString filename="test.root", Bool_t isDebug = 0){
  
  gSystem->Load("libbeam.so");

  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "Raw ROOTFile " << filename << " is opened. " << endl;
  Ssiz_t first_t = filename.Last('/') +1; // if a slash is not there, return 0.
  Ssiz_t last_t = filename.Last('.');
  Int_t length_t = last_t - first_t;
  TString prefix_t = filename(first_t,length_t);
  //__________________________________________________________________________________
  TString output_filename = prefix_t + "_gem.root";
  TFile* rf_rec = TFile::Open(output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;
  // Create Trees 
  TTree* tree_rec = new TTree("Rec","Rec"); // Reconstructed Tree
  //And Build  Reconstruction Branch
  vector <double> vCharge_x, vCharge_y;
  vector <double> vPosition_x, vPosition_y;
  vector <double> vWidth_x, vWidth_y;
  Int_t nHits;
  tree_rec->Branch("charge_x",&vCharge_x);
  tree_rec->Branch("charge_y",&vCharge_y);
  tree_rec->Branch("position_x",&vPosition_x);
  tree_rec->Branch("position_y",&vPosition_y);
  tree_rec->Branch("width_x",&vWidth_x);
  tree_rec->Branch("width_y",&vWidth_y);

  tree_rec->Branch("nHits",&nHits);
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
	// // effectively zero suppression
	if(bgStrip->GetADCsum()>THRESHOLD*sqrt(6)*rms[iProj][ich]){ 
	  bgStrip->Process();
	  bgProjection[iProj]->AddStrip(bgStrip);
	}
      } // End channel loop
      bgProjection[iProj]->Process();
      // bgProjection[iProj]->PlotResults(prefix_t,ievt);
    } // End Projection Loop

    BeamGEMPlane *bgPlane1 = new BeamGEMPlane("gem1");
    bgPlane1->AddProjectionX(bgProjection[0]);
    bgPlane1->AddProjectionY(bgProjection[1]);
    // bgPlane2->AddProjectionX(bgProjection[2]);
    // bgPlane2->AddProjectionY(bgProjection[3]);
    bgPlane1->Process();
    // bgPlane2->Process(); 

    vCharge_x = bgPlane1->GetChargeX();
    vCharge_y = bgPlane1->GetChargeY();
    vPosition_x = bgPlane1->GetPositionX();
    vPosition_y = bgPlane1->GetPositionY();

    vWidth_x = bgPlane1->GetWidthX();
    vWidth_y = bgPlane1->GetWidthY();

    nHits = bgPlane1->GetNHits();
    if(nHits==1){
      bgPlane1->PlotResults(prefix_t,ievt);

    }
    //** Fill this entry
    tree_rec->Fill();
  } // End Event loop

  //__________________________________________________________________________________
  
  rf_raw->Close();
  tree_rec->Write();
  rf_rec->Close();

  return 0;
}
