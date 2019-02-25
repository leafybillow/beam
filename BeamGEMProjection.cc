#include "BeamGEMProjection.h"
#include "BeamGEMStrip.h"
#include "BeamParameters.h"

#include "TCanvas.h"
#include "TText.h"

#include <iostream>

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
  baseline_mean = CalculateMean(vBaseline);
  baseline_rms = CalculateRMS(vBaseline);
  baseline_rms = sqrt( pow(baseline_rms,2) - pow(baseline_mean,2));

  overall_mean = CalculateMean(vStat);
  overall_rms = CalculateRMS(vStat);
  overall_rms = sqrt(pow(overall_rms,2) - pow(overall_mean,2));

  charge_sum =0;
  
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
	aCluster.fCharge = ProcessCharge( vecRange[iCluster]);
	aCluster.fPosition = ProcessCentroid( vecRange[iCluster]);
	aCluster.fWidth = ProcessWidth(vecRange[iCluster]);
	aCluster.peak =  ProcessSplitCheck(vecRange[iCluster]);
	if( (aCluster.peak).size()<=1)
	  aCluster.fSplit = 0;
	else
	  aCluster.fSplit = (aCluster.peak).size()-1;
	
	aCluster.pRange = vecRange[iCluster];
	  
	vClusters.push_back( aCluster);
	charge_sum += aCluster.fCharge;
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
  vHits.clear();
  std::vector< ACluster>::iterator it = vClusters.begin();
  while(it!=vClusters.end()){
    AHit aHit;
    if( (*it).fSplit==0){
      aHit.fPosition = (*it).fPosition;
      aHit.fCharge = (*it).fCharge;
      aHit.fWidth = (*it).fWidth;
      //aHit.fRes = ProcessResolution((*it).pRange);
      vHits.push_back(aHit);
    }
    else{
      int npeaks = (*it).fSplit+1;
      double sum = 0;
      for(int i=0; i<npeaks;i++){
    	int ibin = (*it).peak[i];
    	sum += h_proj->GetBinContent(ibin);
      }
      for(int i=0; i<npeaks;i++){
    	int ibin = (*it).peak[i];
    	double peak_val = h_proj->GetBinContent(ibin);
    	double ratio = peak_val/sum;
    	//FIXME: Now isolate multiple hits according to peak height
    	aHit.fPosition = h_proj->GetBinCenter(ibin);
    	aHit.fCharge = ((*it).fCharge)*ratio;
    	aHit.fWidth = ((*it).fWidth)*ratio;
    	vHits.push_back(aHit);
      }
    }
    it++;
  }
  SortHits();
  return 0 ;
}
vector< pair<int,int> > BeamGEMProjection::SearchClusters(){

  vector< pair<int,int> > vecRange;
  double bin_content;
  bool isLock = 0;
  double threshold =0; // Assuming ZeroSuppression has been done in the analysis script
  
  int low=0, up=0;

  // int edge_cut = 5; // cut out false hits at the edge
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
      else{
	for(int i=low;i<up;i++)
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
  if(vHits.size()>1){
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
  bool zsflag = bgGEMStrip->GetZSStatus();
  
  if(!zsflag){ // if not zero suppressed
    if(strip_id>edge_cut && strip_id <nStrips-edge_cut)
      h_proj->SetBinContent(strip_id,ampl);
  }
  
  if(zsflag)
    vBaseline.push_back(ampl);
  
  h_raw->SetBinContent(strip_id,ampl);
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
  if(vClusters.size()>1){
    std::vector< ACluster>::iterator it = vClusters.begin();
    std::vector< ACluster>::iterator it_next = it+1;
    int isCrossTalk = 0;

    while(it!=vClusters.end()){
      while(it_next!=vClusters.end()){
	isCrossTalk = TestCrossTalk(*it, *it_next);
	if(isCrossTalk){
	  int bin_begin = ((*it_next).pRange).first;
	  int bin_end = ((*it_next).pRange).second;
	  
	  for(int i=bin_begin;i<bin_end;i++)
	    h_proj->SetBinContent(i,0);
	  
	  it_next = vClusters.erase(it_next);
	}
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
  if(j.fCharge > (i.fCharge)*xtalk_threshold)
    isCrossTalk = 0;
  else {
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

vector<int> BeamGEMProjection::ProcessSplitCheck(pair<int,int> prRange){
  
  // Adapted from Jlab-TreeSearch::GEMPlane::Decode()
  int iter = prRange.first;
  double max_val = h_proj->GetBinContent(iter);
  double max_pos = iter;
  vector<int> peak;
  int end = prRange.second;
  int width = prRange.second - prRange.first;

  int kMaxSize = width_threshold;
  double frac = split_frac; 
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
	if(cur_val>max_val){
	  max_val = cur_val;
	  max_pos = iter;
	}
	else if( (max_val-cur_val) >stability || iter == end){
	  eStatus = kFindMin;
	  min_val = cur_val;
	  peak.push_back(max_pos);
	  //	  continue;
	}
	break;
      case kFindMin:
	if(cur_val<min_val)
	  min_val = cur_val;
	else if (cur_val-min_val> stability){
	  eStatus = kFindMax;
	  max_val = cur_val;
	  max_pos = iter;
	  min_counts ++;
	  //	  continue;
	}
	break;
      } // End of switch
    } // End of iter while-loop
  } // End of ... I know
  return peak;  // No splitting found 
}


double BeamGEMProjection::CalculateMean(vector<double> aVector){
  double sum = 0;
  int counts = 0;
  vector<double>::iterator iter = aVector.begin();
  while(iter!=aVector.end()){
    sum += *iter;
    counts++;
    iter++;
  }
  if(counts!=0){
    double mean = sum /counts;
    return mean;
  }
  else
    return 0;
}

double BeamGEMProjection::CalculateRMS(vector<double> aVector){
  double sum = 0;
  int counts = 0;
  vector<double>::iterator iter = aVector.begin();
  while(iter!=aVector.end()){
    double buff = *iter;
    sum += (buff*buff);
    counts++;
    iter++;
  }
  if(counts!=0){
    double rms = sqrt(sum /counts);
    return rms;
  }
  else
    return 0;
}


int BeamGEMProjection::PostProcess(){
  // FIXME: to-do
  return 0;
}
