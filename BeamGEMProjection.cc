#include "BeamGEMProjection.h"
#include "BeamGEMStrip.h"
#include <iostream>
#include "TCanvas.h"
ClassImp(BeamGEMProjection);


BeamGEMProjection::BeamGEMProjection(TString Name, Int_t nch)
  :vBGStrips(NULL),vHits(NULL),nHits(-1),isSplit(0){
  nStrips = nch;
  strProjName = Name;
  Init();
}
BeamGEMProjection::~BeamGEMProjection(){

}

void BeamGEMProjection::Init(){
  double pitch = 0.4 ; // 400 um = 0.4 mm
  double length = pitch*nStrips; // (nStrips-1)+0.5+0.5 = nStrips :)

  h_proj = new TH1D("",Form("Projection on %s",strProjName.Data()),
		    nStrips,-length/2.0,length/2.0);
}

int BeamGEMProjection::Process(){

  if(CheckNStrips()==1){ // if Failed to match nstrip
    return 1;
  }

  vector< pair<int,int> > vecRange=SearchClusters();
  AHit aHit;
  int nClusters = vecRange.size();
  for(int iCluster=0; iCluster <nClusters;iCluster++){
    aHit.fCharge = ProcessCharge( vecRange[iCluster]);
    aHit.fPosition = ProcessCentroid( vecRange[iCluster]);
    aHit.fMpl = ProcessMultiplicity(vecRange[iCluster]);
    aHit.fRes = ProcessResolution(vecRange[iCluster]);

    vHits.push_back(aHit);
  }

  nHits = nClusters; // FIXME: just for now, we will check splitting 
  SortHits();

  return 0;
}

vector< pair<int,int> > BeamGEMProjection::SearchClusters(){

  vector< pair<int,int> > vecRange;
  double bin_content;
  int mpl_cut = 2;
  bool isLock = 0;
  double threshold =0; // Zero suppression has been done in the main script
  
  int low=0, up=0;
  for(int iStrip=0; iStrip<nStrips; iStrip++){
    bin_content = h_proj->GetBinContent(iStrip+1);
    if(bin_content>threshold && isLock==0){
      isLock = 1;
      low= iStrip+1;
    }
    if((bin_content<=threshold||iStrip==nStrips) && isLock==1){
      isLock = 0;
      up = iStrip+1;
      if((up-low+1)>mpl_cut){
	vecRange.push_back( make_pair(low,up) );
      }
    }
  }
  return vecRange;
}

 // Compute total charge and its centroid 
double BeamGEMProjection::ProcessCentroid( pair<int,int> prRange){
  int low = prRange.first;
  int up = prRange.second;
  
  double total_charge = h_proj->Integral(low,up);
  double moment =0;  // sum q*x
  double q,x;
  for(int ibin=low;ibin<up;ibin++){
    q = h_proj->GetBinContent(ibin);
    x = h_proj->GetBinCenter(ibin);
    moment+=q*x;
  }
  double pos =  moment/total_charge;
  return pos;
}

double BeamGEMProjection::ProcessCharge( pair<int,int> prRange){
  int low = prRange.first;
  int up = prRange.second;
  double total_charge = h_proj->Integral(low,up);
  return total_charge;
}

int BeamGEMProjection::ProcessMultiplicity(pair<int,int> prRange){
  int mpl = prRange.second - prRange.first +1;
  return mpl;
}

double BeamGEMProjection::ProcessResolution(pair<int,int> prRange){
  double res = 115; // 400/sqrt(12) ,unit:  um

  int mpl = prRange.second - prRange.first +1;
  //This is a crude approach to spatial resolution
  res = res/ mpl; 
  return res;
}

void BeamGEMProjection::SortHits(){
  // Sort Hits by Charge Amplitude.
  // We do this to prepare for correlation matching in Plane level 
  // insert sort is used here
  AHit aHit_buff;
  for(int iHit=1; iHit<nHits; iHit++){
    aHit_buff = vHits[iHit];

    int jHit = iHit-1;
    while(jHit>=0 && vHits[jHit].fCharge > aHit_buff.fCharge){
      vHits[jHit+1]=vHits[jHit];
      jHit = jHit-1;
    }
    vHits[jHit+1]=aHit_buff;
  }
  // idiot check
  // for(int iHit=0; iHit<nHits-1; iHit++){
  //   if( vHits[iHit].fCharge > vHits[iHit+1].fCharge)
  //     std::cout << "Sorting went wrong! " << std::endl;
  // }
}

int BeamGEMProjection::CheckNStrips(){
  //FIXME
  // char cProjName[] = strProjName[0];
  // if((strcmp(cProjName,"x")== 0 && nStrips == 256)
  //    ||(strcmp(cProjName,"y")==0 && nStrips==512)){

  //   return 0;
  // }
  // else{
  //   std::cout << "Error: Number of Strips mismatches the Projection belonged to" << std::endl;
  //   return 1;
  // }
  return 0;
}

void BeamGEMProjection::AddStrip(BeamGEMStrip* bgGEMStrip){

  double ampl = bgGEMStrip->GetAmplitude();
  int strip_id = bgGEMStrip->GetStripID();
  h_proj->SetBinContent(strip_id,ampl);

}

void BeamGEMProjection::PlotResults(TString runName, int ievt){
  TCanvas *c1 = new TCanvas("","c1", 800,400);
  c1->cd();
  h_proj->Draw();
  h_proj->GetYaxis()->SetTitle("Charge(ADC counts)");
  h_proj->GetXaxis()->SetTitle("Position(mm)");
  h_proj->SetTitle( Form("Projection %s ,  %d Cluster(s) ",strProjName.Data(),nHits));
  c1->SaveAs( Form("%s-%s-evt%d.pdf",runName.Data(),strProjName.Data(), ievt) );

  delete c1;
}
