void PrintCounts(){


  TFile *rootfile = TFile::Open("rootfiles/run_411_analyzed.root");
  TTree *tree = rootfile->Get("Rec");

  Int_t gem1_counts;
  Int_t gem2_counts;

  TBranch* branch_qdc = tree->GetBranch("qdc");
  TLeaf* us_hi  =branch_qdc->GetLeaf("us_hi");
  
  tree->SetBranchAddress("gem1.nHits",&gem1_counts);
  tree->SetBranchAddress("gem2.nHits",&gem2_counts);

  // Int_t nentries = 200; //tree->GetEntries();

  // for(Int_t ievt =0; ievt< nentries; ievt++){
  //   tree->GetEntry(ievt);
  //   Int_t qdc_raw =  us_hi->GetValue();
  //   Double_t qdc_count = (qdc_raw-1300)/140.0;
    
  //   cout<< ievt << "\t"
  // 	<< qdc_raw<< "\t"
  // 	<< qdc_count << "\t"
  // 	<< gem1_counts << "\t"
  // 	<< gem2_counts << "\t"
  // 	<< endl;
  // }
  
  Int_t nentries = tree->GetEntries();
  double *evt = new double[nentries];
  for(Int_t ievt =0; ievt< nentries; ievt++){
    evt[ievt] = ievt;
  }
  
  // Int_t npt = tree->Draw("TMath::Nint((qdc.us_hi-1300)/140) -gem2.nHits","","goff");
  Int_t npt = tree->Draw("gem1.nHits -gem2.nHits","","goff");

  TGraph *g1 = new TGraph(npt,evt,tree->GetV1());
  g1->SetMarkerStyle(20);
  g1->Draw("ALP");
  g1->SetTitle("(GEM 1 count - GEM 2 count) vs Event number");
  g1->GetXaxis()->SetTitle("event number");
}
