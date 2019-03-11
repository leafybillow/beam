#include <iostream>
#include "TMath.h"
#include "BeamGEMPlane.h"
#include "BeamGEMProjection.h"


ClassImp(BeamGEMPlane);

BeamGEMPlane::BeamGEMPlane()
  :fPos_x(),fPos_y(),
   fCharge_x(),fCharge_y(),
   fWidth_x(),fWidth_y(),
   fSplit_x(),fSplit_y(),
   fCorelation(),
   vHitsMask_x(),vHitsMask_y(),
   nHits(-1),
   bgProjX(NULL),bgProjY(NULL),
   my_id(-1)
{}

BeamGEMPlane::~BeamGEMPlane(){}

int BeamGEMPlane::Process(){
  if(CheckProjections()==0){
    nHits = Reconstruct();
    if(nHits!=0)
      CollectResults();
    return 0;
  }
  else
    return 1; // Fail PlaneName Check
}

int BeamGEMPlane::Reconstruct(){
  // some assumption is made to make it works. -TY
  int nHits_x = bgProjX->GetNHits();
  int nHits_y = bgProjY->GetNHits();

  if(nHits_x==0 || nHits_y == 0)
    return 0;
  else if (nHits_x==nHits_y){
    correlator aCorrelator;
    aCorrelator.xHits = xHits;
    aCorrelator.yHits = yHits;
    vCorrelator.push_back(aCorrelator);
    return nHits_x;
  }
  else if(nHits_x> nHits_y){
    int key_y ;
    int skip_key = 0;
    int diff = nHits_x- nHits_y+1;
    for(int iy=0;iy<nHits_y;iy++){
      key_y = ( 1 << iy);
      vector<int> vdummy;
      vector<int> vKey_x = GenerateKeys(nHits_x,diff,0,vdummy) ;
      vector<int>::iterator it_keyx = vKey_x.begin();
      correlator candidate = GenerateCorrelator(*it_keyx,key_y);
      int key_candidate = *it_keyx;
      while( it_keyx != vKey_x.end() ){
	if( ((*it_keyx)&skip_key)>0){
	  it_keyx++;
	  continue;
	}
	correlator aCorrelator = GenerateCorrelator(*it_keyx, key_y);
	if(aCorrelator.charge_distance < candidate.charge_distance){
	  cout << aCorrelator.charge_distance << endl;
	  candidate = aCorrelator;
	  key_candidate = *it_keyx;
	}
	it_keyx ++;
      }
      skip_key |= key_candidate;
      diff -= ((candidate.xHits).size()-(candidate.yHits).size());

      UpdateCorrelator(candidate);
      vCorrelator.push_back(candidate);
    }
    int n_rec = 0;
    vector<correlator>::iterator itc=vCorrelator.begin();
    while(itc!=vCorrelator.end()){
      n_rec += ((*itc).xHits).size();
      itc++;
    }
    cout <<n_rec << endl;
    return n_rec;
  }
  else if(nHits_x< nHits_y){
    int key_x ;
    int skip_key = 0;
    int diff = nHits_y- nHits_x+1;
    for(int ix=0;ix<nHits_x;ix++){
      key_x = ( 1 << ix);
      vector<int> vdummy;
      vector<int> vKey_y = GenerateKeys(nHits_y,diff,0,vdummy) ;
      vector<int>::iterator it_keyy = vKey_y.begin();
      correlator candidate = GenerateCorrelator(key_x,*it_keyy);
      int key_candidate = *it_keyy;
      while( it_keyy != vKey_y.end() ){
	if( ((*it_keyy)&skip_key)>0){
	  it_keyy++;
	  continue;
	}
	correlator aCorrelator = GenerateCorrelator(key_x, *it_keyy);
	if(aCorrelator.charge_distance < candidate.charge_distance){
	  cout << aCorrelator.charge_distance << endl;
	  candidate = aCorrelator;
	  key_candidate = *it_keyy;
	}
	it_keyy ++;
      }
      skip_key |= key_candidate;
      diff -= ((candidate.yHits).size()-(candidate.xHits).size());

      UpdateCorrelator(candidate);
      vCorrelator.push_back(candidate);
    }
    int n_rec = 0;
    vector<correlator>::iterator itc=vCorrelator.begin();
    while(itc!=vCorrelator.end()){
      n_rec += ((*itc).yHits).size();
      itc++;
    }
    cout <<n_rec << endl;
    return n_rec;
  }
  else
    return 0;
}

vector<int> BeamGEMPlane::GenerateKeys(int nlength,int nOccupied,int skip,
				       vector<int> prev_keys){
  if(nOccupied <=0)
    return prev_keys;
  
  vector<int> new_keys;
  if(prev_keys.size()==0){
    int idigi =0;
    while(idigi<nlength){
      new_keys.push_back(1<<idigi);
      idigi++;
    } 
  }
  else{
    new_keys = prev_keys;
    vector<int>::iterator iter_prev = prev_keys.begin();
    iter_prev += skip;
    while(iter_prev!=prev_keys.end()){
      int nshift = 0;
      while ( ((*iter_prev>>nshift) & 1)!=1)
	nshift ++;
      int idigi =0;
      while(idigi<nshift){
	new_keys.push_back(*iter_prev|(1<<idigi));
	idigi++;

      }
      iter_prev++;
    }
    skip+= prev_keys.size();    
  }
  return GenerateKeys(nlength-1,nOccupied-1,skip,new_keys);
}

correlator BeamGEMPlane::GenerateCorrelator(int key_x, int key_y){
  int nbits_x = floor(TMath::Log2(key_x))+1;
  int nbits_y = floor(TMath::Log2(key_y))+1;

  correlator aCorrelator;
  vector< AHit> xHits_buff;
  vector< AHit> yHits_buff;

  for(int ibit =0 ;ibit<nbits_x ;ibit++){
    int key =(key_x>>ibit)&1;
    if(key)
      xHits_buff.push_back(xHits[ibit]);
  }

  for(int ibit =0 ;ibit<nbits_y ;ibit++){
    int key =(key_y>>ibit)&1;
    if(key)
      yHits_buff.push_back(yHits[ibit]);
  }
  
  aCorrelator.xHits = xHits_buff;
  aCorrelator.yHits = yHits_buff;
  
  EvalChargeDistance(aCorrelator);
  return aCorrelator;
}


void BeamGEMPlane::UpdateCorrelator( correlator &aCorrelator){

  int nhits_x = aCorrelator.xHits.size();
  int nhits_y = aCorrelator.yHits.size();

  if(nhits_y==1 && nhits_x>1){
    
    double total_charge_x = aCorrelator.charge_sum_x;
    double total_charge_y = aCorrelator.charge_sum_y;
    
    double position_y = aCorrelator.yHits[0].fPosition;
    double width_y = aCorrelator.yHits[0].fWidth;
    aCorrelator.yHits.clear();
    // re-distribute hits
    for(int i=0; i<nhits_x; i++){
      double ratio = aCorrelator.xHits[i].fCharge/ total_charge_x;
      AHit aHit;
      aHit.fPosition = position_y;
      aHit.fWidth = width_y;
      aHit.fCharge = ratio*total_charge_y;
      aCorrelator.yHits.push_back(aHit);
    }
    
  }
  else if(nhits_x==1 && nhits_y>1){
    double total_charge_x = aCorrelator.charge_sum_x;
    double total_charge_y = aCorrelator.charge_sum_y;
    
    double position_x = aCorrelator.xHits[0].fPosition;
    double width_x = aCorrelator.xHits[0].fWidth;
    aCorrelator.xHits.clear();
    // re-distribute hits
    for(int i=0; i<nhits_y; i++){
      double ratio = aCorrelator.yHits[i].fCharge/ total_charge_y;
      AHit aHit;
      aHit.fPosition = position_x;
      aHit.fWidth = width_x;
      aHit.fCharge = ratio*total_charge_x;
      aCorrelator.xHits.push_back(aHit);
    }
  }
  
}

void BeamGEMPlane::EvalChargeDistance( correlator &aCorrelator){

  double xcharge =0;
  double ycharge =0;

  vector<AHit>::iterator itx = aCorrelator.xHits.begin();
  vector<AHit>::iterator ity = aCorrelator.yHits.begin();

  while(itx!= aCorrelator.xHits.end()){
    xcharge += (*itx).fCharge;
    itx++;
  }

  while(ity!= aCorrelator.yHits.end()){
    ycharge += (*ity).fCharge;
    ity++;
  }
  aCorrelator.charge_sum_x = xcharge;
  aCorrelator.charge_sum_y = ycharge;
  aCorrelator.charge_distance = sqrt(pow(xcharge-ycharge,2));
}


void BeamGEMPlane::CollectResults(){

  vector<correlator>::iterator it = vCorrelator.begin();
  while(it!=vCorrelator.end()){

    vector<AHit>::iterator itx = (*it).xHits.begin();
    vector<AHit>::iterator ity = (*it).yHits.begin();
    
    while(itx!=(*it).xHits.end()){
      fCharge_x.push_back( (*itx).fCharge);
      fPos_x.push_back( (*itx).fPosition);
      fWidth_x.push_back( (*itx).fWidth);
      itx++;
    }

    while(ity!=(*it).yHits.end()){
      fCharge_y.push_back( (*ity).fCharge);
      fPos_y.push_back( (*ity).fPosition);
      fWidth_y.push_back( (*ity).fWidth);
      ity++;
    }

    it++;
  }
}

int BeamGEMPlane::CheckProjections(){
  TString strName_Y = bgProjY->GetProjName();
  TString strName_X = bgProjX->GetProjName();

  if(strName_Y.Contains("y") && strName_X.Contains("x"))
    return 0;
  else{
    std::cerr << "Error:"
	      <<__FUNCTION__
	      <<" Incorrect Projections pair"
	      << std::endl;
    return 1;
  }
}

void BeamGEMPlane::AddProjectionX(BeamGEMProjection* bgProj){
  bgProjX = bgProj;
  xHits = bgProjX->GetHits();
}

void BeamGEMPlane::AddProjectionY(BeamGEMProjection* bgProj){
  bgProjY = bgProj;
  yHits = bgProjY->GetHits();
}

void BeamGEMPlane::PlotResults(TString runName, int ievt){

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(1,2);

  c1->cd(1);
  bgProjY->GetRawHist()->Draw();

  TVirtualPad *c2 = c1->cd(2);
  c2->Divide(2,1);
  c2->cd(1);
  bgProjX->GetRawHist()->Draw();
  
  c1->SaveAs(Form("%s-%s-evt-%d.png",
		  runName.Data(),
		  strPlaneName.Data(),
		  ievt));
  
  delete c1;
}

