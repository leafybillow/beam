void Chain(){
  
  TChain* chain = new TChain("Rec");
  chain->Add("rootfiles/run_132_analyzed.root");
  chain->Add("rootfiles/run_136_analyzed.root");
  chain->Add("rootfiles/run_138_analyzed.root");
  chain->Add("rootfiles/run_142_analyzed.root");
  chain->Add("rootfiles/run_143_analyzed.root");
  chain->Add("rootfiles/run_144_analyzed.root");
  chain->Add("rootfiles/run_145_analyzed.root");
  
  TCanvas *c1 =new TCanvas("c1","c1",800,800);

  c1->cd();
  chain->Draw("det1_qdc_hr>>h0(150,800,1500)","isNoisy==0");
  chain->Draw("det1_qdc_hr>>h1(150,800,1500)","nPrimaries==1 && isNoisy==0","goff");
  chain->Draw("det1_qdc_hr>>h2(150,800,1500)","nPrimaries==2&& isNoisy==0","goff");
  chain->Draw("det1_qdc_hr>>h3(150,800,1500)","nPrimaries==3 && isNoisy==0","goff");

  // chain->Draw("det1_qdc_hr>>h1(200,800,1800)","nTracks==1","goff");
  // chain->Draw("det1_qdc_hr>>h2(200,800,1800)","nTracks==2","goff");
  // chain->Draw("det1_qdc_hr>>h3(200,800,1800)","nTracks==3","goff");
  
  h0->SetLineStyle(2);
  h0->SetLineWidth(2);
  h1->SetLineColor(2);  h1->SetLineWidth(2);
  h2->SetLineColor(4);  h2->SetLineWidth(2);
  h3->SetLineColor(6);  h3->SetLineWidth(2);
  h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
    
  TCanvas *c2 =new TCanvas("c2","c2",800,800);
  c2->cd();
  chain->Draw("(det1_qdc_hr-875)/5.78>>hfit(120,-2,100)","nTracks==1 && isNoisy==0","goff");
  // chain->Draw("det1_qdc_hr-875>>hfit(150,-10,600)","nTracks==1 && isNoisy==0","goff");

  TF1 *flangau = new TF1("flangau",langaufun,0,400,4);

  double par[4];
  flangau->SetParName(0,"Sigma_Landau");
  flangau->SetParName(1,"MPV");
  flangau->SetParName(2,"Integral");
  flangau->SetParName(3,"Sigma_Gaussian");

  // par[0] = 0.1;
  // par[1] = 140;
  // par[2] = hfit->Integral(1,400); // per bin = 1
  // par[3] = 30;

  par[0] = 0.1;
  par[1] = 25;
  par[2] = hfit->Integral(1,80); // per bin = 1
  par[3] = sqrt(par[1]);

  flangau->SetParameters(par);
  hfit->Draw();
  hfit->Fit(flangau,"R+","",10,60);

}


Double_t langaufun(Double_t *x, Double_t *par) {
   //Fit parameters:
   //par[0]=Width (sca1le) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = 0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;

      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
      return (par[2] * step * sum * invsq2pi / par[3]);
}
























