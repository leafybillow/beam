void EvalCorrelation(Int_t run_number){

  TFile* rootfile = TFile::Open(Form("./rootfiles/run_%d_analyzed.root",run_number));
  TTree* rec = (TTree*)rootfile->Get("Rec");
  
  rec->SetAlias("width"," (det2_qdc_hr-1275)/80.0");  
  rec->SetAlias("QDC"," (det1_qdc_hr-1300)");
  rec->SetAlias("correlation","(det-1.662)*(nPrimaries-1.308)/(1.035*1.345)");
  
  // rec->SetAlias("det"," (det1_qdc_hr-875)/130.0");
  // rec->SetAlias("correlation","(det-1.477)*(nPrimaries-1.428)/(0.9423*0.9077)");
  TString cut = "Entry$>=0 && Entry$<200";
  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(1,2);
  c1->cd(1);
  rec->Draw("nPrimaries:Entry$",cut,"*l");
  c1->cd(2);
  rec->Draw("QDC:Entry$",cut,"l*");

  // rec->Draw("correlation:nPrimaries","det2_qdc_hr>1500","*");

  // rec->Draw("correlation:det","det2_qdc_hr>1500 && nPrimaries==0","*");
  // rec->Draw("correlation:det","det2_qdc_hr>1500 && nPrimaries==1","* same");
  // rec->Draw("correlation:det","det2_qdc_hr>1500 && nPrimaries==2","* same");
  // rec->Draw("correlation:det","det2_qdc_hr>1500 && nPrimaries==3","* same");

  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  c2->cd();
  rec->Draw("QDC:nPrimaries>>h2(10,-0.5,9.5,100,0,1000)","","COLZ ");
  h2->SetTitle("QDC Signal v.s. GEM nPrimaries");
  h2->GetYaxis()->SetTitle(" QDC Signal");
  h2->GetXaxis()->SetTitle(" Number of Primary Tracks");
  // rec->Draw("det:width>>(10,-0.5,50,100,0,10)","nPrimaries==0","COLZ ");

  // rec->Draw("det1_y:det1_x","nPrimaries==1","COLZ");
  // rec->Draw("det:nPrimaries>>(10,-0.5,9.5,100,0,10)","det2_qdc_hr>1500 && correlation<1 && correlation>0 ","COLZ");
  // rec->Draw("(det-1.477)/0.9423:(nPrimaries-1.428)/0.9077","","prof");

}
