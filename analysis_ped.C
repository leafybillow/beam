 // This is a script to generate a rootfile containing pedestals 
 // Last Update Jan. 2019  TY

#define N_GEM 2; // Number of GEMS

Int_t analysis_ped(TString filename="test.root", Bool_t isDebug = 0){

  TFile* rf_raw = TFile::Open(filename.Data());
  TTree* tree_raw = (TTree*)rf_raw->Get("T");
  cout << "ROOTFile " << filename << " is opened. " << endl;

  TString output_filename = filename.ReplaceAll(".root","_ped.root");
  TFile* rf_ped = TFile::Open(output_filename,"RECREATE");
  cout << "ROOTFile " << output_filename << " is recreated. " << endl;


  // * Pre-Analysis GEM : Get Pedestal From Gaussian Fit
  const int ngem = N_GEM;
  const int napv = 6;  // number of apvs for each GEM
  const int nch = 128; // number of channel on each APV cards
  TH1D* hped_mean[ngem][napv];
  TH1D* hped_rms[ngem][napv];

  for(int igem=0; igem<ngem;igem++){
    for(int iapv=0; iapv<napv;iapv++){
      hped_mean[igem][iapv] = new TH1D(Form("hped_mean_%d_%d",igem,iapv),
				       Form("Pedestal Mean, GEM %d, APV %d",igem, iapv),
				       128,-0.5,127.5);
      hped_rms[igem][iapv] = new TH1D(Form("hped_rms_%d_%d",igem,iapv),
				      Form("Pedestal RMS, GEM %d, APV %d",igem, iapv),
				      128,-0.5,127.5);
    }
  }

  TString apv_tag[napv] = {"x","x",
			   "y","y","y","y"};
  Int_t base[napv] = {0, 128, 
		      0, 128 , 256, 384};
  Double_t ped_mean; //averaged by 6
  Double_t ped_rms; // averaged by sqrt(6)
  
  TString draw_text, hist_name;
  TH1D* h_fit = new TH1D("h_fit","Buffer historgram for pedestal fit",1e4,-0.5,2e4-0.5);

  cout << "Calculating Pedestals... Be Patient" << endl;
  for(int igem =0; igem< ngem; igem++){
    for(int iapv=0; iapv< napv; iapv++){
      for(int ich=0; ich<nch;ich++){
	draw_text = Form("sbs.gems.%s%d.adc_sum[%d]",
			 apv_tag[iapv].Data(),
			 igem+1,
			 ich+base[iapv]); // stripe number
	tree_raw->Draw(Form("%s >> h_fit", draw_text.Data()),"","goff");
	PedestalFit(h_fit, ped_mean, ped_rms);
	
	hped_mean[igem][iapv]->SetBinContent(ich+1,ped_mean);
	hped_rms[igem][iapv]->SetBinContent(ich+1,ped_rms);
      } 
    } 
  }
  
  rf_raw->Close();

  //Write  objects  
  for(int igem =0; igem<ngem;igem++){
    for(int iapv=0; iapv<napv; iapv++){
      hped_mean[igem][iapv]->Write();
      hped_rms[igem][iapv]->Write();
    }
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
  h_ped->Fit("f_gaus","QN","",bincenter-2*rms,bincenter+2*rms);
  
  mean = f_gaus->GetParameter(1)/6.0; // averaged by 6
  sigma  = f_gaus->GetParameter(2)/sqrt(6);  // averaged by sqrt(6), assuming samples are not correlated

 }
