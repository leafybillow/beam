#include "BeamGEMProjection.h"
#include "BeamGEMStrip.h"
#include <iostream>
#include "TCanvas.h"
ClassImp(BeamGEMProjection);
#define WIDTH_CUT 1 ; // Rejecting single active channel 

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

  if(CheckNStrips()==1) // if it fails to match nstrip
    return 1;

  else{

    if(h_proj->GetEntries()>0){
      vector< pair<int,int> > vecRange=SearchClusters();
      AHit aHit;
      int nClusters = vecRange.size();
      for(int iCluster=0; iCluster <nClusters;iCluster++){

	aHit.fCharge = ProcessCharge( vecRange[iCluster]);
	aHit.fPosition = ProcessCentroid( vecRange[iCluster]);
	aHit.fWidth = ProcessWidth(vecRange[iCluster]);
	aHit.fRes = ProcessResolution(vecRange[iCluster]);

	vHits.push_back(aHit);
      }

      nHits = nClusters; // FIXME: just for now, we will check splitting 
      SortHits();
    } // Pass if h_proj has non-zero entries
    else 
      nHits = 0;

    h_proj->SetTitle( Form("Projection %s ,  %d Cluster(s) Found ",strProjName.Data(),nHits));

    return 0;
  } // Pass CheckNStrips()

}

vector< pair<int,int> > BeamGEMProjection::SearchClusters(){

  vector< pair<int,int> > vecRange;
  double bin_content;
  int width_cut = WIDTH_CUT;
  bool isLock = 0;
  double threshold =0; // Assuming ZeroSuppression has been done in the analysis script
  
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
      if((up-low)>width_cut){
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
  double moment =0;  // sum of q*x
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

int BeamGEMProjection::ProcessWidth(pair<int,int> prRange){
  int width = prRange.second - prRange.first;
  return width;
}

double BeamGEMProjection::ProcessResolution(pair<int,int> prRange){
  double res = 115; // 400/sqrt(12) ,unit:  um

  int width = prRange.second - prRange.first +1;
  //This is a crude approach to spatial resolution
  res = res/ width; 
  return res;
}

void BeamGEMProjection::SortHits(){
  // Sort Hits by Charge Amplitude.
  // We do this to prepare for correlation matching in Plane level 
  // insertion sort is used here
  AHit aHit_buff;
  for(int iHit=1; iHit<nHits; iHit++){
    aHit_buff = vHits[iHit];
    int jHit = iHit-1;
    while(jHit>=0 && vHits[jHit].fCharge < aHit_buff.fCharge){
      vHits[jHit+1]=vHits[jHit];
      jHit = jHit-1;
    }
    vHits[jHit+1]=aHit_buff;
  }
  // idiot check
  for(int iHit=0; iHit<nHits-1; iHit++){
    if( vHits[iHit].fCharge < vHits[iHit+1].fCharge)
      std::cout << "Sorting went wrong! " << std::endl;
  }
}

int BeamGEMProjection::CheckNStrips(){
  
  if( (strProjName.Contains("x") && nStrips == 256)
      ||(strProjName.Contains("y") && nStrips== 512)){
    return 0;
  }
  else{
    std::cout <<"Line:" <<  __LINE__ << " " 
	      <<__FUNCTION__ << " "
	      << "Error: Number of Strips mismatches the Projection belonged to"
	      << std::endl;
    std::cout << "strProjName: " << strProjName << std::endl;
    std::cout << "nStrips: " << nStrips << std::endl;
    return 1;
  }
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

  c1->SaveAs( Form("%s-%s-evt%d.pdf",runName.Data(),strProjName.Data(), ievt) );

  delete c1;
}


int BeamGEMProjection::TestCrossTalk(int iHit1, int iHit2){

  // int strip[4] = {32,88,118,127} ;  unit  number strips
  double target[4] = {12.8, 35.2, 47.2, 50.8}; // target = strip*0.4 unit: mm
  
  double position1 = vHits[iHit1].fPosition;
  double width = vHits[iHit1].fWidth * 0.5; // just need half the width
  double position2 = vHits[iHit2].fPosition;

  int myapv1  = position1 /50.8;
  int myapv2  = position2 /50.8; // returns an integer representing apv id

  double separation = fabs(position2-position1);
  bool isCrossTalk = 0; 
  int i=0;

  if( myapv1 != myapv2 ){
    // APV ids mismatch with each other : 
    // two positions are crossing the border of 2 APVs , so it is not cross-talk
    return 0;
  }
  else{
    
    while(isCrossTalk==0 && i!=4 ){

      if( fabs(separation-target[i])<=width ) 
	isCrossTalk = 1;
      
      i++;
    }

    if(isCrossTalk)
      return 1;
    else
      return 0;
  }

}

void BeamGEMProjection::UpdateHits( vector< int> vHitsMask){

  // 1. Remove cross talk or noisy clusters from vHits array
  // 2. Update Histogram title, indicating number of hits identified
  
  int nMask = vHitsMask.size();
  int sumMask = 0;

  if(nMask!= nHits){
    std::cerr << "Error : " 
	      << __FILE__ << ":"
	      << __FUNCTION__ <<  "()::"
	      << " Mismatched number of masks and hits"
	      << std::endl;
  }
  else{

    std::vector< AHit>::iterator iHits = vHits.begin();
    for(int iMask=0; iMask<nMask; iMask++){
      
      if(vHitsMask[iMask]==1){ // A Good Hit
	sumMask +=1;
	iHits++;
      }
      else
	iHits = vHits.erase(iHits);
    }
    
  }
  
  if(sumMask == vHits.size())
    nHits = vHits.size();
  else{
    std::cerr << "Error : " 
	      << __FILE__ << ":"
	      << __FUNCTION__ <<  "()::"
	      << " Mismatched masks sum  and  number of passed  hits"
	      << std::endl;
  }
  
  TString title = h_proj->GetTitle();
  TString append = Form(", %d hit(s) confirmed", nHits);
  h_proj->SetTitle(title+append);

}
