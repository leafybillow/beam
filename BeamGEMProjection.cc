#include "BeamGEMProjection.h"
#include <iostream>
ClassImp(BeamGEMProjection);


BeamGEMProjection::BeamGEMProjection(TH1D* h, double *rms, TString Name)
  :vHits(NULL),nStrips(-1),nHits(-1),isSplit(0){
  h_proj =h;
  fRMS = rms;
  strName = Name;
  Init();
}
BeamGEMProjection::~BeamGEMProjection(){

}

void BeamGEMProjection::Init(){
  nStrips = h_proj->GetNbinsX();
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
  int mpl_cut = 3;
  bool isLock = 0;
  
  // bool isAccept = 0;
  int low=0, up=0;
  for(int iStrip=0; iStrip<nStrips; iStrip++){
    bin_content = h_proj->GetBinContent(iStrip+1);
    if(bin_content>0 && isLock==0){
      isLock = 1;
      low=bin_content;
    }
    if((bin_content<=0||iStrip==nStrips) && isLock==1){
      isLock = 0;
      up = bin_content;
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
  // Bubble sort is used here
  AHit aHit_buff;
  for(int iHit=1; iHit<nHits; iHit++){
    for(int jHit=0; jHit<iHit; jHit++){
      if(vHits[jHit].fCharge > vHits[iHit].fCharge){
	aHit_buff = vHits[iHit];
	vHits[iHit] = vHits[jHit];
	vHits[jHit] = aHit_buff;
      }
    }
  }
}

int BeamGEMProjection::CheckNStrips(){
  if((strName=="x" && nStrips == 256) || (strName=="y" &&nStrips==512)){
    return 0;
  }
  else{
    std::cout << "Error: Number of Strips mismatches the Projection belonged to" << std::endl;
    return 1;
  }

}
