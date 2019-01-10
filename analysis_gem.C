#define N_GEM 2; // Number of GEMS
#define THRESHOLD 3.0; // threshold for signal, e.g. 3 means 3*sigma 
#define N_CH_CMN 100;  // Number of strips used for common mode noise calculation
#define PITCH 400.0;   // pitch width, unit: um
#define CUT_MPL 1;     // cut on multiplicity, reject single wire hits (most likely noises or cross-talks)
#define CUT_COR 1.0;     // cut on correlation, based on correlation studies

Int_t analysis_gem(TString filename="test.root", Bool_t isDebug = 0){
  
  gSystem->Load("libbeam.so");

  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "Raw ROOTFile " << filename << " is opened. " << endl;
  //__________________________________________________________________________________
  TString output_filename = filename.ReplaceAll(".root","_analyzed.root");
  TFile* rf_rec = TFile::Open(output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;
  // Create Trees 
  TTree* tree_rec = new TTree("Rec","Rec"); // Reconstructed Tree
  //And Build  Reconstruction Branch

  //__________________________________________________________________________________
  // GEM Configuration Parameters
  const Int_t ngem = N_GEM;
  const Int_t napv = 6;  // number of apvs for each GEM
  const Int_t nch = 128; // number of channel on each APV cards

  TString apv_tag[napv] = {"x","x",
			   "y","y","y","y"};
  Int_t base[napv] = {0, 128, 
		      0, 128 , 256, 384};

  Double_t cmn[ngem][napv];
  Int_t nch_cmn = N_CH_CMN; // number of channels used for computing  common mode

  for(Int_t igem =0;igem<ngem;igem++){
    for(Int_t iapv=0;iapv<napv;iapv++){
      cmn[igem][iapv] = 0.0 ;
    }
  }

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
    Double_t* id_strip; // channel-to-strip mapping array
    Double_t nChannel;  // Note: nch in SBS-offline decoder is a float number

    // identification for projeciton
    Int_t id_Proj; // starts from 0 to 3; {"x1","y1","x2","y2"}
  };
  struct StcProjection Proj[4];

  for(Int_t iProj=0;iProj<nProj;iProj++){

    Proj[iProj].id_strip = new Double_t[ sizeArray[iProj] ];
    Proj[iProj].adc_sum = new Double_t[ sizeArray[iProj] ];
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

    tree_raw->SetBranchAddress(Form("sbs.gems.%s.adc_sum",strProj[iProj].Data()),
			       Proj[iProj].adc_sum);

    tree_raw->SetBranchAddress(Form("sbs.gems.%s.nch",strProj[iProj].Data()),
			       &(Proj[iProj].nChannel));
  }

  
  //__________________________________________________________________________________
  //*Event loop, reconstruction
  Int_t nentries = 10;//tree_raw->GetEntries();
  for(Int_t ievt=0;ievt<nentries;ievt++){
    tree_raw->GetEntry(ievt);
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
	for(Int_t ich=0; ich<nChannel;ich++){
	  
	  if( Proj[iProj].adc_sum[ich] > 200){
	    Double_t arADC[6];  	  // *** Transpose ADC Samples
	    for(Int_t iadc=0 ;iadc<nadc;iadc++){
	      arADC[iadc] = Proj[iProj].adc[iadc][ich];
	    } // End sample loop
	    
	    BeamGEMStrip *bgStrip = new BeamGEMStrip(arADC);
	    bgStrip->Process();
	    std::cout << bgStrip->GetRawAmplitude()
	    	      << "\t"
	    	      << bgStrip->GetFitAmplitude() 
	    	      <<"\t" 
	    	      << bgStrip->GetADCsum() 
	    	      <<"\t" 
	    	      << bgStrip->GetADCsum()/ bgStrip->GetRawAmplitude()
	    	      <<"\t" 
	    	      << bgStrip->GetFitAmplitude()/bgStrip->GetRawAmplitude()<<std::endl;
	    
	  } // if this channel survives Zero Suppression

	} // End channel loop
      } // End number of channel check
    } // End Projection Loop

    //** Fill this entry
    //tree_rec->Fill();
    //** Option: output plots
  }  // End Event loop

  //__________________________________________________________________________________
  rf_raw->Close();
  tree_rec->Write();
  rf_rec->Close();
  return 0;
}

// User Define Functions
Double_t EvalCommonMode(Int_t *arr_integral){
  
}

Double_t CommonModeCorrection(Double_t *samples, Double_t fCMN){
  
}
