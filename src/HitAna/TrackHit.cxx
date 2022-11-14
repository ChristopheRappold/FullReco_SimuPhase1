#include "ConstantParameter.hh"
#include "TrackHit.hh"
#include "FiberHitXUV.hh"
#include "TString.h"
#include "FiberHitAna.hh"
#include "FiberTrackAna.hh"
#include "MDCHitAna.hh"
#include "ParaManager.hh"
#include <iostream>
#include <math.h>

TrackHit::TrackHit(){}
TrackHit::~TrackHit(){}

void TrackHit::Print(){
  for( auto v: _fiber_cont) v->Print();
  for( auto v: _mdc_cont)   v->Print();
}


void TrackHit::Add(FiberHitAna *hit_fiber){
  _fiber_cont.emplace_back(hit_fiber);
}

void TrackHit::Add(MDCHitAna *hit_mdc){
  _mdc_cont.emplace_back(hit_mdc);
}

void TrackHit::DeleteDupMDC(){

  for(int i=0; i<(int)_mdc_cont.size(); ++i){
    int lay_i = _mdc_cont[i]->GetLay();
    int tot_i = _mdc_cont[i]->GetTot();
    double dl_i = _mdc_cont[i]->GetDl();
    double dif_i = _mdc_cont[i]->GetDif();
    for(int j=(int)_mdc_cont.size()-1; j>=0; --j){
      int lay_j = _mdc_cont[j]->GetLay();
      int tot_j = _mdc_cont[j]->GetTot();
      double dl_j = _mdc_cont[j]->GetDl();
      double dif_j = _mdc_cont[j]->GetDif();
      //if(lay_i==lay_j && tot_j < tot_i){
      //if(lay_i==lay_j && dl_i < dl_j){
      if(lay_i==lay_j && fabs(dif_i) < fabs(dif_j)){
        _mdc_cont.erase(_mdc_cont.begin() + j);
        if(j<i) i--;
      }
    }
  }

}




FiberHitAna*  TrackHit::GetFiberHit(int i){
  if(i<0 || i >= (int)_fiber_cont.size()){
    std::cout << "i is wrong" << std::endl;
    return 0;
  }
  return _fiber_cont[i];
}

MDCHitAna*  TrackHit::GetMDCHit(int i){
  if(i<0 || i >= (int)_mdc_cont.size()){
    std::cout << "i is wrong" << std::endl;
    return 0;
  }
  return _mdc_cont[i];
}

void TrackHit::SetTrack(FiberTrackAna *track, ParaManager *par){
  _x = track->GetXtgt();
  _y = track->GetYtgt();
  _a = track->GetA();
  _b = track->GetB();
  _z = par->fiber_tgt_pos_z;
}
