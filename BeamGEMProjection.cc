#include "BeamGEMProjection.h"
#include "BeamGEMStrip.h"
#include <iostream>
#include "TCanvas.h"
ClassImp(BeamGEMProjection);
#define WIDTH_CUT 1 ; // Rejecting single active channel 
#define THRESHOLD_WIDTH 3 ; // Threshold to examine oversize cluster, unit: number of strip
#define SPLIT_FRAC 0.1 ;
#define PITCH 0.4 ; // unit mm,  = 400 um

BeamGEMProjection::BeamGEMProjection()
  :vBGStrips(),vHits(),vClusters(),
   nHits(-1),nClusters(-1),
   strProjName(),
   h_proj(),h_raw()
{
}
BeamGEMProjection::BeamGEMProjection(TString Name, Int_t nch)
  :vBGStrips(),vHits(),vClusters(),
   nHits(-1),nClusters(-1)
{
  nStrips = nch;
  strProjName = Name;
  Init();
}
BeamGEMProjection::~BeamGEMProjection(){

}

void BeamGEMProjection::Init(){
  double pitch = PITCH ; // 
  double length = pitch*nStrips; // (nStrips-1)+0.5+0.5 = nStrips :)

  h_proj = new TH1D("",Form("Projection on %s",strProjName.Data()),
		    nStrips,-length/2.0,length/2.0);
  
  h_raw = new TH1D("",Form("Projection on %s",strProjName.Data()),
		    nStrips,-length/2.0,length/2.0);
  h_raw->GetYaxis()->SetTitle("Charge(ADC counts)");
  h_raw->GetXaxis()->SetTitle("Position(mm)");
}
int BeamGEMProjection::Process(){

  int status = CoarseProcess();
  if(status==0)
    status = FineProcess();
    
  return status;
}
int BeamGEMProjection::CoarseProcess(){

  if(CheckNStrips()==1) // if it fails to match nstrip
    return 1;//Failed
  else{
    if(h_proj->GetEntries()>0){
      vector< pair<int,int> > vecRange = SearchClusters();
      ACluster aCluster;
      nClusters = vecRange.size();
      
      for(int iCluster=0; iCluster <nClusters;iCluster++){
	  aCluster.fCharge = ProcessCharge( vecRange[iCluster]);
	  aCluster.fPosition = ProcessCentroid( vecRange[iCluster]);
	  aCluster.fWidth = ProcessWidth(vecRange[iCluster]);
	  aCluster.fSplit =  ProcessSplitCheck(vecRange[iCluster]);
	  aCluster.pRange = vecRange[iCluster];
	  
	  vClusters.push_back( aCluster);
      }
      SortClusters();
      RejectCrossTalk();
      // update number of clusters after cross talk rejection
      nClusters = vClusters.size();
    }// Pass if h_proj has non-zero entries
    else
      nClusters = 0;

    // summarize nHits counting splitting peaks
    nHits = nClusters;
    for(int iCluster =0; iCluster<nClusters;iCluster++){
      nHits += vClusters[iCluster].fSplit;
    }
    
    if(nClusters==0)
      h_raw->SetTitle( Form("Projection %s ,  %d Cluster(s) Found",
			     strProjName.Data(),
			     nClusters));
    else if(nClusters>0)
      h_raw->SetTitle( Form("Projection %s ,  %d Cluster(s) Found, %d Hit(s) Identified",
			     strProjName.Data(),
			     nClusters,
			     nHits));
    return 0; // OK
  } // Pass CheckNStrips()
}
int BeamGEMProjection::FineProcess(){
  //FIXME: TO-DO
  if(nHits==1){ // Right now only process single non-splitting cluster hits
    AHit aHit;
    aHit.fPosition = vClusters[0].fPosition;
    aHit.fCharge = vClusters[0].fCharge;
    aHit.fWidth = vClusters[0].fWidth;
    aHit.fRes = ProcessResolution(vClusters[0].pRange);
    
    vHits.push_back(aHit);
  }
  return 0 ;
}
vector< pair<int,int> > BeamGEMProjection::SearchClusters(){

  vector< pair<int,int> > vecRange;
  double bin_content;
  int width_cut = WIDTH_CUT;
  bool isLock = 0;
  double threshold =0; // Assuming ZeroSuppression has been done in the analysis script
  
  int low=0, up=0;

  int edge_cut = 5; // cut out false hits at the edge
  int start = edge_cut;
  int end = nStrips-edge_cut;
  
  for(int iStrip= start; iStrip<end; iStrip++){
    bin_content = h_proj->GetBinContent(iStrip+1);
    if(bin_content>threshold && isLock==0){
      isLock = 1;
      low= iStrip+1;
    }
    if((bin_content<=threshold||iStrip==end) && isLock==1){
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
void BeamGEMProjection::SortClusters(){
  // Sort Cluster by Charge Amplitude.
  // We do this to prepare for cross talk check 
  // insertion sort is used here
  ACluster aCluster_buff;
  for(int i=1; i<nClusters; i++){
    aCluster_buff = vClusters[i];
    int j = i-1;
    while(j>=0 && vClusters[j].fCharge < aCluster_buff.fCharge){
      vClusters[j+1]=vClusters[j];
      j = j-1;
    }
    vClusters[j+1]=aCluster_buff;
  }
  // idiot check
  for(int i=0; i<nClusters-1; i++){
    if( vClusters[i].fCharge < vClusters[i+1].fCharge)
      std::cout << "Sorting went wrong! " << std::endl;
  }
}

void BeamGEMProjection::SortHits(){
  // Sort Hits by Charge Amplitude.
  // We do this to prepare for correlation matching in Plane level 
  // insertion sort is used here
  AHit aHit_buff;
  for(int i=1; i<nHits; i++){
    aHit_buff = vHits[i];
    int j = i-1;
    while(j>=0 && vHits[j].fCharge < aHit_buff.fCharge){
      vHits[j+1]=vHits[j];
      j = j-1;
    }
    vHits[j+1]=aHit_buff;
  }
  // idiot check
  for(int i=0; i<nHits-1; i++){
    if( vHits[i].fCharge < vHits[i+1].fCharge)
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
}

void BeamGEMProjection::AddStrip(BeamGEMStrip* bgGEMStrip){

  double ampl = bgGEMStrip->GetAmplitude();
  int strip_id = bgGEMStrip->GetStripID();
  bool zsflag = bgGEMStrip->GetZSStatus();
  
  if(!zsflag){ // if not zero suppressed
    h_proj->SetBinContent(strip_id,ampl);
  }

  h_raw->SetBinContent(strip_id,ampl);
}

void BeamGEMProjection::PlotResults(TString runName, int ievt){
  TCanvas *c1 = new TCanvas("","c1", 800,400);
  c1->cd();
  h_raw->Draw();

  c1->SaveAs( Form("%s-%s-evt%d.pdf",runName.Data(),strProjName.Data(), ievt) );

  delete c1;
}

void BeamGEMProjection::RejectCrossTalk(){
  if(vClusters.size()>1){
    std::vector< ACluster>::iterator it = vClusters.begin();
    std::vector< ACluster>::iterator it_next = it+1;
    int isCrossTalk = 0;

    while(it!=vClusters.end()){
      while(it_next!=vClusters.end()){
	isCrossTalk = TestCrossTalk(*it, *it_next);
	if(isCrossTalk)
	  it_next = vClusters.erase(it_next);
	else
	  it_next++;
      }
      it++;
      it_next = it+1;
    } // end of outer loop
  } // if only one or less, nothing need to be done.
}

int BeamGEMProjection::TestCrossTalk(ACluster i, ACluster j){

  // Not a very good idea for oversize cluster
  // It is tedious but more safe.

  int isCrossTalk = 0;
  // induced cluster is usually relatively small
  // a center value is good enough

  int strip2_center = (j.pRange.second +j.pRange.first)/2.0;
  int myapv2 = floor(strip2_center/128); // APV #(0-3)
  strip2_center = strip2_center%128; // reduced to APV strip #(0-127)

  int strip2_neighbor_lo = lower_neighbor[strip2_center] + myapv2*128;
  int strip2_neighbor_hi = higher_neighbor[strip2_center]+ myapv2*128;

  // the range of main cluster
  int strip1_lo = i.pRange.first;
  int strip1_up = i.pRange.second;
  // Calculate the APV id separately,
  // because the main cluster could cross the border of two APVs
  int myapv1_lo = floor(strip1_lo/128);
  int myapv1_up = floor(strip1_up/128);
  
  // Since induced cluster sits in a relative small range,
  // it would be easier to start searching from the induced one 
  if( (myapv2-myapv1_up)*(myapv2-myapv1_lo) == 0) {
    // at least one is equal and then trigger cross talk test
    int a = (strip2_neighbor_hi - strip1_lo)*(strip2_neighbor_hi - strip1_up);
    int b = (strip2_neighbor_lo - strip1_lo)*(strip2_neighbor_lo - strip1_up);
    if ( a<0 || b<0)
      isCrossTalk =1;
    else
      isCrossTalk =0;
  }
  else
    isCrossTalk = 0; // it is not

  return isCrossTalk;
}

// int BeamGEMProjection::TestCrossTalk_v1(int iHit1, int iHit2){
// //A quick way to test cross talk
//   // int strip[4] = {32,88,118,127} ;  unit  number strips
//   double target[4] = {12.8, 35.2, 47.2, 50.8}; // target = strip*0.4 unit: mm
  
//   double position1 = vHits[iHit1].fPosition;
//   double width = vHits[iHit1].fWidth * 0.5; // just need half the width
//   double position2 = vHits[iHit2].fPosition;

//   int myapv1  = ceil(position1 /50.8);
//   int myapv2  = ceil(position2 /50.8); // returns an integer representing apv id

//   double separation = fabs(position2-position1);
//   bool isCrossTalk = 0; 
//   int i=0;

//   if( myapv1 != myapv2 ){
//     // APV ids mismatch with each other : 
//     // two positions are crossing the border of 2 APVs , so it is not cross-talk
//     return 0;
//   }
//   else{
//     while(isCrossTalk==0 && i!=4 ){
//       if( fabs(separation-target[i])<=width ){
// 	isCrossTalk = 1;
// 	break;
//       }
//       i++;
//     }
//     if(isCrossTalk)
//       return 1;
//     else
//       return 0;
//   }
// }

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
  
  if(sumMask == (int)vHits.size())
    nHits = vHits.size();
  else{
    std::cerr << "Error : " 
	      << __FILE__ << ":"
	      << __FUNCTION__ <<  "()::"
	      << " Mismatched masks sum  and  number of passed  hits"
	      << std::endl;
  }
  
  TString title = h_raw->GetTitle();
  TString append = Form(", %d hit(s) confirmed", nHits);
  h_raw->SetTitle(title+append);

}

int BeamGEMProjection::ProcessSplitCheck(pair<int,int> prRange){
  
  // Adapted from Jlab-TreeSearch::GEMPlane::Decode()
  int iter = prRange.first;
  double max_val = h_proj->GetBinContent(iter);
  int end = prRange.second;
  int width = prRange.second - prRange.first;

  int kMaxSize = THRESHOLD_WIDTH;
  double frac = SPLIT_FRAC;
  double frac_up = 1.0 + frac;
  double frac_down = 1.0 -frac;

  double cur_val ; //current bin content  
  double min_val;
  int min_counts = 0; 
  enum EStatus {kFindMax=1,kFindMin};
  EStatus eStatus = kFindMax;
  if( width >kMaxSize ){
    while(iter!=end){
      iter++;
      cur_val = h_proj->GetBinContent(iter);
      switch(eStatus){
      case kFindMax:
	if(cur_val>max_val)
	  max_val = cur_val;
	else if( cur_val < max_val*frac_down){
	  eStatus = kFindMin;
	  min_val = cur_val;
	  continue;
	}
	break;
      case kFindMin:
	if(cur_val<min_val)
	  min_val = cur_val;
	else if (cur_val> min_val*frac_up){
	  eStatus = kFindMax;
	  max_val = cur_val;
	  min_counts ++;
	  continue;
	}
	break;
      } // End of switch
    } // End of iter while-loop
    return min_counts; // Number of local minimum found
  } // End of ... I know
  else
    return 0;  // No splitting found 
}
