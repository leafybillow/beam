void GoldTracks(){
  TString rf_path = getenv("BEAM_ROOTFILES");
  TFile* rootfile = TFile::Open(rf_path+"run_419_analyzed.root");

  TCanvas *c1 =new TCanvas("c1","c1",600,600);
  c1->cd();
  Rec->Draw("det1_qdc_hr","!isNoisy && isGoldenTrack");
  TH1D* h_buff = (TH1D*)c1->FindObject("htemp");
  h_buff->SetLineColor(kRed);


  TCanvas *c2 =new TCanvas("c2","c2",600,600);
  c2->cd();
  Rec->Draw("det1_qdc_hr:det1_x","!isNoisy && isGoldenTrack","prof");


  TCanvas *c4 =new TCanvas("c4","c4",500,800);
  c4->cd();
  c4->SetLogz();
  Rec->Draw("det1_y:det1_x>>(100,-50,50,200,-100,100)","(det1_qdc_hr-1280)*(isGoldenTrack && !isNoisy)","COLZ");

}







