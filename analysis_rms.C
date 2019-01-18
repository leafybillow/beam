// This is a script to generate a header file and a root file 
// containing RMSs after common mode correction
// Last Update Jan. 2019  TY

Int_t analysis_rms(TString filename="test.root", Bool_t isDebug = 0){

  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "ROOTFile " << filename << " is opened. " << endl;

  Ssiz_t first_t = filename.Last('/') +1; // if a slash is not there, return 0.
  Ssiz_t last_t = filename.Last('.');
  Int_t length_t = last_t - first_t;
  TString prefix_t = filename(first_t,length_t);
  TString output_filename = prefix_t + "_rms.root";
  TFile* rf_rms = TFile::Open("rootfiles/"+output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;

  // Header for GEM
  TString header_filename = Form("DBfiles/analysis_rms.h_%s",prefix_t.Data());
  FILE *header_file = fopen(header_filename.Data(),"w");
  
  // Insert a comment line 
  fprintf(header_file,"\n");
  fprintf(header_file," //// RMS arrays %s", prefix_t.Data());
  
  const int nproj = 4;
  TH1D* hped_mean[nproj];
  TH1D* hped_rms[nproj];

  TString strProj[nproj]={"x1","y1","x2","y2"};
  Int_t sizeArray[nproj]={256, 512, 256, 512}; 


  for(int iproj=0; iproj<nproj;iproj++){
    int nbins =sizeArray[iproj];
    hped_mean[iproj] = new TH1D(Form("hped_mean_%d",iproj),
				Form("Pedestal Mean vs strip,  %s",strProj[iproj].Data()),
				nbins,-0.5,nbins-0.5);
    hped_rms[iproj] = new TH1D(Form("hped_rms_%d",iproj),
			       Form("Pedestal RMS vs strip,  %s",strProj[iproj].Data()),
			       nbins,-0.5,nbins-0.5);
  }


  // Retrieve Channel Mapping from rootfiles ........
  double gem1_xstrip_id[256];
  double gem1_ystrip_id[512];
  double gem2_xstrip_id[256];
  double gem2_ystrip_id[512];
  double* strip_id[nproj] = { gem1_xstrip_id, gem1_ystrip_id,
			      gem2_xstrip_id, gem2_ystrip_id }; 
  // strip_id type should be an integer, but ....it was output as a float
  // Hardcoded by hand ,may need to think of a better solution to do this.
  // Anyway, I need a pointer to load strip map
  tree_raw->SetBranchAddress("sbs.gems.x1.strip",gem1_xstrip_id);
  tree_raw->SetBranchAddress("sbs.gems.y1.strip",gem1_ystrip_id);
  tree_raw->SetBranchAddress("sbs.gems.x2.strip",gem2_xstrip_id);
  tree_raw->SetBranchAddress("sbs.gems.y2.strip",gem2_ystrip_id);
  
  tree_raw->GetEntry(1); // in order to load strip map to these array
  // Caution: And this needs to be fixed if ZeroSuppression is on in SBS-offline
  
  Double_t ped_mean; //averaged by 6
  Double_t ped_rms; // averaged by sqrt(6)
  
  TString draw_text, hist_name;
  TH1D* h_fit = new TH1D("h_fit","Buffer historgram for pedestal fit",1e3,-5e3,5e3);

  cout << "Calculating rms... Be patient" << endl;

  for(int iproj=0; iproj<4; iproj++){

    fprintf(header_file,"\n");
    fprintf(header_file,"double rms_%s[%d]={",strProj[iproj].Data(),sizeArray[iproj]);

    int nch = sizeArray[iproj];
    for(int ich=0; ich<nch;ich++){
      draw_text = Form("sbs.gems.%s.adc_sum[%d]-6*sbs.gems.%s.common_mode[%d]",
		       strProj[iproj].Data(),ich,
		       strProj[iproj].Data(),ich);
      tree_raw->Draw(Form("%s >> h_fit", draw_text.Data()),"","goff");
      PedestalFit(h_fit, ped_mean, ped_rms);
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


  printf("Done ! \n");
  rf_raw->Close();

  fclose(header_file);
  //Write  objects  
  for(int iproj=0;iproj<nproj;iproj++){
      hped_mean[iproj]->Write();
      hped_rms[iproj]->Write();
  }
  rf_rms->Close();

  return 0;
}

// User Define Functions
void PedestalFit(TH1D *h_ped, Double_t &mean, Double_t &sigma){

  Int_t bin_max = h_ped->GetMaximumBin();
  Double_t bincenter = h_ped->GetBinCenter(bin_max);
  Double_t bin_content_max = h_ped->GetBinContent(bin_max);
  Double_t rms = 50.0; // An initial guess

  Double_t par[3]; 
  par[0] = bin_content_max;
  par[1] = bincenter;
  par[2] = rms; 

  TF1 *f_gaus = new TF1("f_gaus","gaus",0,5e4);
  f_gaus->SetParameters(par);
  h_ped->Fit("f_gaus","QNR","",bincenter-rms,bincenter+rms);

  mean = f_gaus->GetParameter(1);
  sigma  = f_gaus->GetParameter(2);
  h_ped->Fit("f_gaus","QNR","",mean-sigma,mean+sigma);

  mean = f_gaus->GetParameter(1)/6.0; // averaged by 6
  sigma  = f_gaus->GetParameter(2)/sqrt(6);  // averaged by sqrt(6), assuming samples are not correlated
}
