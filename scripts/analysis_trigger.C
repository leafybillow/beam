#include "analysis_rms.h"
#define THRESHOLD 3.0; // threshold for signal, e.g. 3 means 3*sigma 

Int_t analysis_trigger(TString filename="test.root", Bool_t isDebug = 0){
  
  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "Raw ROOTFile " << filename << " is opened. " << endl;
  Ssiz_t first_t = filename.Last('/') +1; // if a slash is not there, return 0.
  Ssiz_t last_t = filename.Last('.');
  Int_t length_t = last_t - first_t;
  TString prefix_t = filename(first_t,length_t);
  //__________________________________________________________________________________
  TString output_filename = prefix_t + "_triggertime.root";
  TFile* rf_rec = TFile::Open("rootfiles/"+output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;
  // Create Trees 
  TTree* tree_rec = new TTree("Rec","Rec"); // Reconstructed Tree
  //And Build  Reconstruction Branch
  Double_t ev_num = 0;
  Double_t coarse_time1,coarse_time2;
  Double_t fine_time;
  tree_rec->Branch("ev_num",&ev_num);
  tree_rec->Branch("coarse_time1",&coarse_time1);
  tree_rec->Branch("coarse_time2",&coarse_time2);
  tree_rec->Branch("fine_time",&fine_time);

  
  tree_raw->SetBranchAddress("sbs.gems.x1.coarse_time1",&coarse_time1);
  tree_raw->SetBranchAddress("sbs.gems.x1.coarse_time2",&coarse_time2);
  tree_raw->SetBranchAddress("sbs.gems.x1.fine_time",&fine_time);
  
  //__________________________________________________________________________________

  //*Event loop, reconstruction
  
  Int_t nentries = tree_raw->GetEntries();
  for(Int_t ievt=0;ievt<nentries;ievt++){
    if(ievt%200==0)
      cout << ievt << " events analyzed"  << endl;

    tree_raw->GetEntry(ievt);
    ev_num ++;

    tree_rec->Fill();
  } // End Event loop

  //__________________________________________________________________________________
  
  rf_raw->Close();
  tree_rec->Write();
  rf_rec->Close();

  return 0;
}
