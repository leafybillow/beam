 // This is a script to generate a rootfile containing pedestals 
 // Last Update Jan. 2019  TY

#define N_GEM 2; // Number of GEMS

Int_t analysis_ped(TString filename="test_raw.root", Bool_t isDebug = 0){

  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "ROOTFile " << filename << " is opened. " << endl;

  TString output_filename = filename.ReplaceAll(".root","_ped.root");
  TFile* rf_ped = TFile::Open(output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;
  // Database for GEM
  TString db_filename = "db_sbs.gems.dat_test";
  gSystem->Exec("cp db_sbs.gems.dat_noped "+db_filename);
  FILE *db_file = fopen(db_filename.Data(),"a");
  gSystem->Exec(Form("ln -sf %s db_sbs.gems.dat",db_filename.Data()));
  
  // Insert a comment line 
  fprintf(db_file,"\n");
  fprintf(db_file,"# Pedestal Database");
  
  // * Pre-Analysis GEM : Get Pedestal From Gaussian Fit

  const int nproj = 4;
  TH1D* hped_mean[nproj];
  TH1D* hped_rms[nproj];

  TString strProj[nproj]={"x1","y1","x2","y2"};
  Int_t sizeArray[nproj]={256, 512, 256, 512}; 

  for(int iproj=0; iproj<nproj;iproj++){
      hped_mean[iproj] = new TH1D(Form("hped_mean_%d",iproj),
				  Form("Pedestal Mean vs strip, projection %s",strProj[iproj].Data()),
				  sizeArray[iproj],-0.5,sizeArray[iproj]-0.5);
      hped_rms[iproj] = new TH1D(Form("hped_rms_%d",iproj),
				 Form("Pedestal RMS vs strip, projection %s",strProj[iproj].Data()),
				 sizeArray[iproj],-0.5,sizeArray[iproj]-0.5);
  }

  Double_t ped_mean; //averaged by 6
  Double_t ped_rms; // averaged by sqrt(6)
  
  TString draw_text, hist_name;
  TH1D* h_fit = new TH1D("h_fit","Buffer historgram for pedestal fit",1e4,-0.5,2e4-0.5);
  
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

  cout << "Calculating Pedestals... Be patient" << endl;
  int counts = 0;
  for(int iproj =0; iproj< nproj; iproj++){

    fprintf(db_file,"\nsbs.gems.%s.ped =",strProj[iproj].Data());
    int nch = sizeArray[iproj];

    for(int ich=0; ich<nch;ich++){
      draw_text = Form("sbs.gems.%s.adc_sum[%d]",
		       strProj[iproj].Data(),
		       ich); // channel number
      tree_raw->Draw(Form("%s >> h_fit", draw_text.Data()),"","goff");
      PedestalFit(h_fit, ped_mean, ped_rms);
	
      // Get Strip id from array
      int strip = strip_id[iproj][ich]; 
      // File-Print pedestals
      if(ich%8==0)
	fprintf(db_file, " \\ \n");

      fprintf(db_file,"%d %.2f ",(int)strip, ped_mean);
      // counts++;
      // if(counts%5==0){
      // 	printf("\r Running %.1f %%  ",counts/total_counts*100);
      // }
      hped_mean[iproj]->SetBinContent(strip+1,ped_mean);
      hped_rms[iproj]->SetBinContent(strip+1,ped_rms);
    } 
  }
  printf("Done ! \n");
  rf_raw->Close();

  fclose(db_file);
  //Write  objects  
  
  for(int iproj=0; iproj<nproj; iproj++){
    hped_mean[iproj]->Write();
    hped_rms[iproj]->Write();
  }

  rf_ped->Close();

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
  h_ped->Fit("f_gaus","QN","",bincenter-rms,bincenter+rms);
  
  mean = f_gaus->GetParameter(1)/6.0; // averaged by 6
  sigma  = f_gaus->GetParameter(2)/sqrt(6);  // averaged by sqrt(6), assuming samples are not correlated

 }
