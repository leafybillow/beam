void MatchKey(Int_t run_number = 412){
  TString rf_path = TString(getenv("BEAM_ROOTFILES"));
  TFile *rootfile = TFile::Open(rf_path + Form("run_%d_analyzed.root",run_number));

  TTree* Rec = (TTree*)rootfile->Get("Rec"); // Reconstructed Tree

  vector<Int_t> GEM_key;
  vector<Int_t> QDC_key;
  Int_t EventLength = 600;
  Int_t KeyLength = EventLength / 2 ;
  
  Int_t nPrimaries ; // Number of Primary Tracks from GEMs
  Double_t qdc_value ; // QDC 
  Double_t qdc_pedestal = 1300 ;
  Double_t threshold = 10 ; // if QDC is less than <threshold> above pedestal 
  Int_t key_true =1 ;
  Int_t key_false =0 ;
  
  Rec->SetBranchAddress("nPrimaries", &nPrimaries);
  Rec->SetBranchAddress("det1_qdc_hr", &qdc_value);
  // Collecting Pattern from GEM and QDC
  for(int ievt =0 ;ievt<KeyLength;ievt++){
    Rec->GetEntry(ievt);
    if(nPrimaries == 0)
      GEM_key.push_back(key_true); // Form Key pattern based on 'zero' event
    else
      GEM_key.push_back(key_false);
  }

  for(int ievt =0 ;ievt<EventLength;ievt++){
    Rec->GetEntry(ievt);
    if(qdc_value-qdc_pedestal<=threshold)
      QDC_key.push_back(key_true);
    else
      QDC_key.push_back(key_false);
  }
  
  // I used GEM key as a reference
  // Now loop over event shifts from 0 to MaxShift
  // and find the best fit between QDC key to GEM key.

  vector<int>::iterator itr_qdc = QDC_key.begin();
  
  Int_t BestGoodness = 0; // Goodness of Key Fit
  Int_t EventShift = 0;
  BestGoodness = GetGoodness(GEM_key, itr_qdc+EventShift); // Initial Guess

  for(int ishift = 0; ishift< KeyLength;ishift++){
    Int_t myGoodness = GetGoodness(GEM_key, itr_qdc+ishift);

    if(myGoodness>BestGoodness){
      BestGoodness=myGoodness;
      EventShift = ishift;
    }
  }
  
  cout << "Event Offset : " << EventShift << endl;
  cout << "Goodness of Fit : " << BestGoodness << endl;
  
}

Int_t GetGoodness(vector<int> key_ref, vector<int>::iterator itr_candidate){

  Int_t ret_Goodness = 0; // Goodness to return

  vector<int>::iterator itr_ref =key_ref.begin();
  
  while(itr_ref!=key_ref.end()){
    ret_Goodness +=( (*itr_ref)*(*itr_candidate) );
    itr_ref++;
    itr_candidate++;
  }

     
  return ret_Goodness;
}
