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

TrackHit::TrackHit(){_mdchit_dif.resize(17); }
TrackHit::~TrackHit(){}



void TrackHit::AddFiber(int layer, int fiberID){
  _fiberhit[layer]=fiberID;
}

void TrackHit::AddPSB(int psbID){
  _psbhit=psbID;
}

void TrackHit::AddPSFE(int psfeID){
  _psfehit=psfeID;
}


void TrackHit::DeleteDupMDC(){
  for(int i=0; i<(int)17; ++i)
  {
    double min_dif=9999.;
    int min_hit=-1;

    for(auto x : _mdchit_dif[i])
    {
      if(x.second < min_dif)
      {
        min_dif=x.second;
        min_hit=x.first;
      }
    }

    _mdchit_bestdif[i]=min_hit;
  }
}




bool TrackHit::IsInclusive(TrackHit* track){

  bool flag = false;

  if( GetNum() <= track->GetNum()) flag = false;
  else{
    bool flag_tmp = true;
    std::vector<int> fiberhit1 = GetFiberHit();
    std::vector<int> fiberhit2 = track->GetFiberHit();
    for(size_t i = 0; i < 6; ++i)
      {
        if(fiberhit2[i] == -1) continue;
        if(fiberhit1[i] != fiberhit2[i]) flag_tmp = false;
      }

    std::vector<int> mdchit1 = GetMDCHit();
    std::vector<int> mdchit2 = track->GetMDCHit();
    for(size_t i = 0; i < 17; ++i)
      {
        if(mdchit2[i] == -1) continue;
        if(mdchit1[i] != mdchit2[i]) flag_tmp = false;
      }

    if(flag_tmp) flag = true;
  }

  return flag;

}

void TrackHit::SetTrack(FiberTrackAna *track, double tgt_posz){
  _x = track->GetXtgt();
  _y = track->GetYtgt();
  _a = track->GetA();
  _b = track->GetB();
  _z = tgt_posz;
}

/*
void TrackHit::ReplaceMDC(std::map<int, std::pair<TVector3, TVector3> > track_pos_mom, ParaManager *par){

  for(int i=0; i<GetNumMDC(); ++i){
    int did   = _mdc_cont[i]->GetDid();
    int layer = _mdc_cont[i]->GetLay();
    int wire  = _mdc_cont[i]->GetWir();
    if(track_pos_mom.find(did)!=track_pos_mom.end() && _mdc_layer_cont[layer].size()>1){
      TVector3 track_pos = track_pos_mom[did].first;
      TVector3 track_mom = track_pos_mom[did].second;
      double dist_min = 9999;
      MDCHitAna *hit_buf;
      for(auto hit: _mdc_layer_cont[layer]){
        TVector3 mdc_pos = hit->GetWirPos(par, track_pos.z());
        TVector3 mdc_dir = hit->GetWirDir(par);
        double dist = -9999;
        //TVector3 buf_vertex = ana->GetVertexPoint( track_pos, mdc_pos, track_mom, mdc_dir, dist);

        double xi=track_pos.x(), yi=track_pos.y(), xo=mdc_pos.x(), yo=mdc_pos.y();
        double ui=track_mom.x()/track_mom.z(), vi=track_mom.y()/track_mom.z();
        double uo=mdc_dir.x()/mdc_dir.z(), vo=mdc_dir.y()/mdc_dir.z();

        double z=((xi-xo)*(uo-ui)+(yi-yo)*(vo-vi))/
          ((uo-ui)*(uo-ui)+(vo-vi)*(vo-vi));
        double x1=xi+ui*z, y1=yi+vi*z;
        double x2=xo+uo*z, y2=yo+vo*z;
        dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

        if(dist < dist_min){
          hit_buf = hit;
          dist_min = dist;
        }
      }

      _mdc_cont.erase(_mdc_cont.begin() + i);
      _mdc_cont.insert(_mdc_cont.begin() + i, hit_buf);

    }
  }
}
*/