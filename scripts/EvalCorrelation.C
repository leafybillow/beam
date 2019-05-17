void EvalCorrelation(Int_t run_number){

  TString rf_path = TString(getenv("BEAM_ROOTFILES"));
  TFile *rootfile = TFile::Open(rf_path + Form("run_%d_analyzed.root",run_number));

  TTree* rec = (TTree*)rootfile->Get("Rec");
  
  rec->SetAlias("QDC","det1_qdc_hr");

  TString cut = "Entry$>=0 && Entry$<200";
  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(1,2);
  c1->cd(1);
  rec->Draw("nPrimaries:Entry$",cut,"*l");
  c1->cd(2);
  rec->Draw("QDC:Entry$",cut,"l*");

  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  c2->cd();
  rec->Draw("QDC:nPrimaries>>h2(10,-0.5,9.5,100,1200,2300)","","COLZ ");
  h2->SetTitle("QDC Signal v.s. GEM nPrimaries");
  h2->GetYaxis()->SetTitle(" QDC Signal");
  h2->GetXaxis()->SetTitle(" Number of Primary Tracks");

}
