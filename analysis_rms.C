// This is a script to generate a header file and a root file 
// containing RMSs after common mode correction
// Last Update Jan. 2019  TY

Int_t analysis_rms(TString filename="test.root", Bool_t isDebug = 0){

  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "ROOTFile " << filename << " is opened. " << endl;

  TString output_filename = filename.ReplaceAll(".root","_rms.root");
  TFile* rf_rms = TFile::Open(output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;

  // Header for GEM
  TString header_filename = "analysis_rms.h";
  FILE *header_file = fopen(header_filename.Data(),"w");
  
  // Insert a comment line 
  fprintf(header_file,"\n");
  fprintf(header_file," //// RMS arrays");
  
  const int nproj = 4;
  TH1D* hped_mean[nproj];
  TH1D* hped_rms[nproj];

  TString strProj[nproj]={"x1","y1","x2","y2"};
  Int_t sizeArray[nproj]={256, 512, 256, 512}; 


  for(int iproj=0; iproj<nproj;iproj++){
    int nbins =sizeArray[iproj];
    hped_mean[iproj] = new TH1D(Form("hped_mean_%d",iproj),
				Form("Pedestal Mean,  %s",strProj[iproj].Data()),
				nbins,-0.5,nbins-0.5);
    hped_rms[iproj] = new TH1D(Form("hped_rms_%d",iproj),
			       Form("Pedestal RMS,  %s",strProj[iproj].Data()),
			       nbins,-0.5,nbins-0.5);
  }

  
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

      hped_mean[iproj]->SetBinContent(ich+1,ped_mean);
      hped_rms[iproj]->SetBinContent(ich+1,ped_rms);
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
