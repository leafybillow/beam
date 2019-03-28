#include "BeamGEMProjection.h"
#include "BeamGEMStrip.h"
#include "BeamParameters.h"

#include "TCanvas.h"
#include "TText.h"

#include <iostream>
#include <bits/stdc++.h>

ClassImp(BeamGEMProjection);

#define PITCH 0.4 ; // unit mm,  = 400 um

BeamGEMProjection::BeamGEMProjection()
  :vBGStrips(),vHits(),vClusters(),
   charge_sum(0),nStrips(0),
   nHits(-1),nClusters(-1),
   strProjName(),
   h_proj(),h_raw(),
   vStat(),vBaseline(),
   baseline_mean(0),baseline_rms(0),overall_mean(0),overall_rms(0)
{
}
BeamGEMProjection::BeamGEMProjection(TString Name, Int_t nch)
  :vBGStrips(),vHits(),vClusters(),
   charge_sum(0),
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

  h_proj->Reset();
  h_raw->Reset();
}

int BeamGEMProjection::Process(){

  int status = CoarseProcess();
  if(status==0)
    status = FineProcess();

  return status;
}
int BeamGEMProjection::CoarseProcess(){
  // Compute baseline RMS and mean;
  sort(vStat.begin(), vStat.end());
  
  baseline_mean = CalculateMean(vStat, 0.1);
  baseline_rms = CalculateRMS(vStat, 0.1);
  baseline_rms = sqrt( pow(baseline_rms,2) - pow(baseline_mean,2));
  
  overall_mean = CalculateMean(vStat,1.0);
  overall_rms = CalculateRMS(vStat,1.0);
  overall_rms = sqrt(pow(overall_rms,2) - pow(overall_mean,2));

  FillProjection();
  
  if(CheckNStrips()==1){ // if it fails to match nstrip
    cerr<< "Failed to match NStrip"<< endl;
    return 1;//Failed
  }
  else{
    if(h_proj->GetEntries()>0){
      vector< pair<int,int> > vecRange = SearchClusters();
      nClusters = vecRange.size();
      for(int iCluster=0; iCluster <nClusters;iCluster++){
	ACluster aCluster;
	aCluster.valley =  FindValleys(vecRange[iCluster]);
	aCluster.fSplit = (aCluster.valley).size();
	aCluster.pRange = vecRange[iCluster];
	  
	vClusters.push_back( aCluster);
      }
      nClusters = vClusters.size();
    }// Pass if h_proj has non-zero entries
    else
      nClusters = 0;

    return 0; // OK
  } // Pass CheckNStrips()
}
int BeamGEMProjection::FineProcess(){
  vHits.clear();
  charge_sum = 0;
  
  std::vector< ACluster>::iterator it = vClusters.begin();
  while(it!=vClusters.end()){
    AHit aHit;
    pair<int,int> pair_range =(*it).pRange;
    pair<int,int> pair_buff;
    vector< pair<int, int> > vec_range ;
    // split one pair to multiplets depending on splitting level
    int nSplit = (*it).fSplit;
    int low = pair_range.first;
    if (nSplit==0)
      vec_range.push_back(pair_range);
    else{
      for(int isplit=0;isplit<nSplit;isplit++){
	pair_buff = make_pair(low,((*it).valley)[isplit]);
	low = ((*it).valley)[isplit];
	vec_range.push_back(pair_buff);
      }
      pair_buff = make_pair(low,pair_range.second);
      vec_range.push_back(pair_buff);
    }
      
    vector< pair<int,int> >::iterator it_pair = vec_range.begin();
    while( it_pair!= vec_range.end() ){
      aHit.fPosition = ProcessCentroid( *it_pair);
      aHit.fHeight = ProcessPeakHeight( *it_pair);
      aHit.fCharge = ProcessCharge( *it_pair);
      aHit.fWidth = ProcessWidth( *it_pair);
      aHit.fRes = ProcessResolution(*it_pair);
      aHit.pRange = *it_pair;
      charge_sum += aHit.fCharge;
      vHits.push_back(aHit);
      it_pair++;
    }
    it++;
  }

  SortHits();
  RejectCrossTalk();
  
  // vector<AHit>::iterator ith = vHits.begin();
  // while(ith!=vHits.end()){
  //   cout << strProjName << ":" ;
  //   cout << (*ith).pRange.first << "-" << (*ith).pRange.second ;
  //   cout << ", charge: ";
  //   cout << (*ith).fCharge << endl;
  //   ith++;
  // }

  // summarize nHits counting splitting peaks

  nHits = vHits.size();	 
  h_raw->SetTitle( Form("Projection %s ,  %d Hit(s) Found",
			strProjName.Data(),
			nHits));

  return 0 ;
}
vector< pair<int,int> > BeamGEMProjection::SearchClusters(){

  vector< pair<int,int> > vecRange;
  double bin_content;
  bool isLock = 0;
  double threshold =0; // Assuming ZeroSuppression has been done in the analysis script
  
  int low=0, up=0;

  int start = edge_cut +1;
  int end = nStrips-edge_cut;
  
  for(int iStrip= start; iStrip<=end; iStrip++){
    bin_content = h_proj->GetBinContent(iStrip);
    
    if(bin_content>threshold && isLock==0){
      isLock = 1;
      low= iStrip;
    }
    if((bin_content<=threshold || iStrip==end) && isLock==1){
      
      if(iStrip!=end){
      	double next = h_proj->GetBinContent(iStrip+1);
      	if(next > threshold){
	  // h_proj->SetBinContent(iStrip, (bin_content+next)*0.5); // some dirty trick
      	  continue;
	}
      }
      
      isLock = 0;
      up = iStrip;
      if((up-low)>width_cut && h_proj->Integral(low,up)>1000){
	  vecRange.push_back( make_pair(low,up) );
      }
      else{
	for(int i=low;i<=up;i++)
	  h_proj->SetBinContent(i,0);
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
  for(int ibin=low;ibin<=up;ibin++){
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
  int width = prRange.second - prRange.first+1;
  return width;
}

double BeamGEMProjection::ProcessPeakHeight(pair<int,int> prRange){
  
  int low = prRange.first;
  int up = prRange.second;
  
  double peak_height = h_proj->GetBinContent(low);
  double bin_val ;
  for(int ibin=low+1;ibin<=up;ibin++){
    bin_val = h_proj->GetBinContent(ibin);
    if(bin_val>peak_height)
      peak_height=bin_val;
  }

  return peak_height;
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
  nHits = vHits.size();
  if(nHits>1){
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
  
  h_raw->SetBinContent(strip_id,ampl);

  if(strip_id>edge_cut && strip_id<(nStrips-edge_cut))
    vStat.push_back(ampl);

}

void BeamGEMProjection::PlotResults(TString runName, int ievt){
  TCanvas *c1 = new TCanvas("","c1", 800,400);
  c1->cd();
  h_raw->Draw();

  TString title = h_raw->GetTitle();
  title = title+Form(", Noise RMS %d",(int)baseline_rms);
  h_raw->SetTitle(title);
  
  TText *text= new TText(0.0,0.95,
			 Form("%s-%s-evt-%d",runName.Data(),strProjName.Data(),ievt));
  text->SetNDC();
  text->Draw("same");

  
  c1->SaveAs( Form("%s-%s-evt%d.png",runName.Data(),strProjName.Data(), ievt) );

  delete c1;
}

void BeamGEMProjection::RejectCrossTalk(){
  if(vHits.size()>1){
    std::vector< AHit>::iterator it = vHits.begin();
    std::vector< AHit>::iterator it_next = it+1;
    int isCrossTalk = 0;

    while(it!=vHits.end()){
      while(it_next!=vHits.end()){
	isCrossTalk = TestCrossTalk(*it, *it_next);
	if(isCrossTalk){
	  ErasePeakFromHist(*it_next);
	  it_next = vHits.erase(it_next);
	}
	else
	  it_next++;
      }
      it++;
      it_next = it+1;
    } // end of outer loop
  } // if only one or less, nothing need to be done.
}

void BeamGEMProjection::ErasePeakFromHist(AHit aHit){
  int bin_low = aHit.pRange.first;
  int bin_up = aHit.pRange.second;
  for(int ibin=bin_low; ibin<=bin_up;ibin++){
    h_proj->SetBinContent(ibin,0);
  }
}

int BeamGEMProjection::TestCrossTalk(AHit i, AHit j){

  // Not a very good idea for oversize cluster
  // It is tedious but more safe.
  int isCrossTalk = 0;
  // induced cluster is usually relatively small
  // a peak position value is good enough
  int begin = j.pRange.first;
  int end = j.pRange.second;
  
  int strip2_center = (begin+end)/2.0-1;
  int myapv2 = floor(strip2_center/128); // APV #(0-3)
  strip2_center = strip2_center%128; // reduced to APV strip #(0-127)

  int strip2_neighbor_lo = lower_neighbor[strip2_center] + myapv2*128;
  int strip2_neighbor_hi = upper_neighbor[strip2_center]+ myapv2*128;

  // the range of main cluster
  int strip1_lo = i.pRange.first-1;
  int strip1_up = i.pRange.second-1;
  // Calculate the APV id separately,
  // because the main cluster could cross the border of two APVs
  int myapv1_lo = floor(strip1_lo/128);
  int myapv1_up = floor(strip1_up/128);
  
  // Since induced cluster sits in a relative small range,
  // it would be easier to start searching from the induced one
  if(j.fHeight > (i.fHeight)*xtalk_threshold)
    isCrossTalk = 0;
  else {
    if( (myapv2-myapv1_up)*(myapv2-myapv1_lo) == 0) {
      // at least one is equal and then trigger cross talk test
      int a = (strip2_neighbor_hi - strip1_lo)*(strip2_neighbor_hi - strip1_up);
      int b = (strip2_neighbor_lo - strip1_lo)*(strip2_neighbor_lo - strip1_up);
      if ( a<=0 || b<=0)
	isCrossTalk =1;
      else
	isCrossTalk =0;
    }
    else
      isCrossTalk = 0; // it is not
  }

  return isCrossTalk;
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

vector<int> BeamGEMProjection::FindValleys(pair<int,int> prRange){
  // Adapted from Jlab-TreeSearch::GEMPlane::Decode()
  int iter = prRange.first;
  double max_val = h_proj->GetBinContent(iter);

  vector<int> valley;
  double valley_pos;
  
  int end = prRange.second;
  
  // int width = prRange.second - prRange.first;
  // int kMaxSize = width_threshold;
  double frac = split_frac; 
  double frac_up = 1.0 + frac;
  double frac_down = 1.0 -frac;

  double cur_val ; //current bin content  
  double min_val;

  int min_counts = 0; 
  enum EStatus {kFindMax=1,kFindMin};
  EStatus eStatus = kFindMax;
  
  while(iter<=end){
    iter++;
    cur_val = h_proj->GetBinContent(iter);

    if(cur_val ==0 && iter!=end)// skip dead channel
      continue;
      
    switch(eStatus){
    case kFindMax:
      if(cur_val>max_val){
	max_val = cur_val;
      }
      else if( cur_val<max_val*frac_down){
	eStatus = kFindMin;
	min_val = cur_val;
	valley_pos = iter;
      }
      break;
    case kFindMin:
      if(cur_val<min_val){
	min_val = cur_val;
	valley_pos = iter;
      }
      else if (cur_val> min_val*frac_up){
	eStatus = kFindMax;
	max_val = cur_val;
	valley.push_back(valley_pos);
	min_counts ++;
      }
      break;
    } // End of switch
  } // End of iter while-loop

  return valley;
}


double BeamGEMProjection::CalculateMean(vector<double> aVector,
					double portion =0.0){
  double sum = 0;
  int counts = 0;

  int nsize = aVector.size();
  int stop = (int)(nsize * (1.0-portion));
  int start = (int)(nsize * portion);

  for(int i=start; i<stop;i++){
    sum += aVector[i];
    counts++;
  }

  if(counts!=0){
    double mean = sum /counts;
    return mean;
  }
  else
    return 0;
}

double BeamGEMProjection::CalculateRMS(vector<double> aVector,
				       double portion = 0.0){
  double sum = 0;
  int counts = 0;
  
  int nsize = aVector.size();
  int stop = (int)(nsize * (1.0-portion));
  int start = (int)(nsize * portion);
  for(int i=start; i<stop;i++){
    sum += aVector[i] * aVector[i];
    counts++;
  }

  if(counts!=0){
    double rms = sqrt(sum /counts);
    return rms;
  }
  else
    return 0;
}

void BeamGEMProjection::FillProjection(){

  for(int i=edge_cut; i<=nStrips-edge_cut; i++){
    double bin_content = h_raw ->GetBinContent(i);
    if( (bin_content - baseline_mean)> zs_threshold*baseline_rms )
      h_proj->SetBinContent(i,bin_content);
  }
  
}

int BeamGEMProjection::PostProcess(){
  // FIXME: to-do
  return 0;
}
