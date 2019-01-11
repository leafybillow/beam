#include "analysis_rms.h"
#define THRESHOLD 3.0; // threshold for signal, e.g. 3 means 3*sigma 

Int_t analysis_gem(TString filename="test.root", Bool_t isDebug = 0){
  
  gSystem->Load("libbeam.so");

  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "Raw ROOTFile " << filename << " is opened. " << endl;
  TString runName = filename.ReplaceAll(".root",""); // I need it .
  //__________________________________________________________________________________
  TString output_filename = runName + "_analyzed.root";
  TFile* rf_rec = TFile::Open(output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;
  // Create Trees 
  TTree* tree_rec = new TTree("Rec","Rec"); // Reconstructed Tree
  //And Build  Reconstruction Branch
  vector <double> vCharge_x, vCharge_y;
  vector <double> vPosition_x, vPosition_y;
  tree_rec->Branch("charge_x",&vCharge_x);
  tree_rec->Branch("charge_y",&vCharge_y);
  tree_rec->Branch("position_x",&vCharge_x);
  tree_rec->Branch("position_y",&vCharge_y);
  //__________________________________________________________________________________
  // GEM Configuration Parameters

  Int_t nProj = 4; // number of projections, 2 GEM *(X+Y) = 4
  Int_t nadc = 6; // number of adc samples

  TString strADC[6]={"adc0","adc1","adc2","adc3","adc4","adc5"};
  TString strProj[4]={"x1","y1","x2","y2"};
  Int_t sizeArray[4]={256, 512, 256, 512}; 
  //__________________________________________________________________________________
  // Initialize EventReader for Raw Tree
  struct StcProjection{
    // containers for data
    Double_t* adc[6];
    Double_t* adc_sum;
    Double_t* common_mode;
    Double_t* id_strip; // channel-to-strip mapping array
    Double_t nChannel;  // Note: nch in SBS-offline decoder is a float number

    // identification for projection
    Int_t id_Proj; // starts from 0 to 3; {"x1","y1","x2","y2"}
  };
  Double_t* rms[4]={rms_x1,rms_y1,rms_x2,rms_y2};
  struct StcProjection Proj[4];

  for(Int_t iProj=0;iProj<nProj;iProj++){

    Proj[iProj].id_strip = new Double_t[ sizeArray[iProj] ];
    Proj[iProj].adc_sum = new Double_t[ sizeArray[iProj] ];
    Proj[iProj].common_mode = new Double_t[ sizeArray[iProj] ];
    Proj[iProj].id_Proj = iProj;

    for(Int_t iadc=0;iadc<nadc;iadc++){
      Proj[iProj].adc[iadc] = new Double_t[ sizeArray[iProj] ];
      TString strBranch = Form("sbs.gems.%s.%s",
			       strProj[iProj].Data(),
			       strADC[iadc].Data());
      tree_raw->SetBranchAddress(strBranch, Proj[iProj].adc[iadc]);
    }
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.strip",strProj[iProj].Data()),
			       Proj[iProj].id_strip);

    tree_raw->SetBranchAddress(Form("sbs.gems.%s.nch",strProj[iProj].Data()),
			       &(Proj[iProj].nChannel));
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.adc_sum",strProj[iProj].Data()),
			       Proj[iProj].adc_sum);
    tree_raw->SetBranchAddress(Form("sbs.gems.%s.common_mode",strProj[iProj].Data()),
			       Proj[iProj].common_mode);


  }
  //__________________________________________________________________________________

  //*Event loop, reconstruction
  Int_t nentries = 50;//tree_raw->GetEntries();
  for(Int_t ievt=0;ievt<nentries;ievt++){
    tree_raw->GetEntry(ievt);

    BeamGEMProjection* bgProjection[4];

    //** Retrieve replayed data from the raw tree
    for(Int_t iProj=0;iProj<nProj;iProj++){
      Int_t nChannel = (Int_t)Proj[iProj].nChannel;
      if ( nChannel != sizeArray[iProj]){       // ***Check number of Channel 
	std::cout << __LINE__ << ":"
		  << "Error: " 
		  << "number of channels  mistached "
		  << " on Projection " << strProj[iProj]<< std::endl;
	std::cout << sizeArray[iProj] << " is expected." << std::endl;
	std::cout << nChannel << " is found." << std::endl;
	break;
      }
      else{  
	bgProjection[iProj] = new BeamGEMProjection( strProj[ Proj[iProj].id_Proj ],
								 (Int_t)Proj[iProj].nChannel);
	for(Int_t ich=0; ich<nChannel;ich++){
	  

	    Double_t arADC[6];  	  // *** Transpose ADC Samples
	    for(Int_t iadc=0 ;iadc<nadc;iadc++){
	      arADC[iadc] = Proj[iProj].adc[iadc][ich] - Proj[iProj].common_mode[ich];
	    } // End sample loop
	    Int_t myStripID = (Int_t)Proj[iProj].id_strip[ich];
	    BeamGEMStrip* bgStrip = new BeamGEMStrip(arADC,myStripID);
	    // Effectively zero suppression
	    if(bgStrip->GetRawAmplitude()>THRESHOLD*rms[iProj][ich]){ 
	      bgStrip->Process();
	      bgProjection[iProj]->AddStrip(bgStrip);
	    }
	} // End channel loop
	bgProjection[iProj]->Process();
	// bgProjection[iProj]->PlotResults(runName,ievt);
      } // End number of channel check
    } // End Projection Loop

    BeamGEMPlane *bgPlane1 = new BeamGEMPlane("gem1");
    BeamGEMPlane *bgPlane2 = new BeamGEMPlane("gem2");
    bgPlane1->AddProjectionX(bgProjection[0]);
    bgPlane1->AddProjectionY(bgProjection[1]);
    bgPlane2->AddProjectionX(bgProjection[2]);
    bgPlane2->AddProjectionY(bgProjection[3]);

    bgPlane1->PlotResults(runName,ievt);
    bgPlane2->PlotResults(runName,ievt);
    //** Fill this entry
    //tree_rec->Fill();
  } // End Event loop

  //__________________________________________________________________________________
  
  rf_raw->Close();
  tree_rec->Write();
  rf_rec->Close();

  return 0;
}
