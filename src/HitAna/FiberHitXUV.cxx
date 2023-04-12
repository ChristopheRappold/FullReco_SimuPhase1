#include "ConstantParameter.hh"
#include "FiberHitAna.hh"
#include "FiberHitXUV.hh"
#include "TString.h"
#include <iostream>
#include <math.h>

FiberHitXUV::FiberHitXUV(double _pos_x, double _pos_y, double _d, FiberHitAna* hit_x, FiberHitAna* hit_u, FiberHitAna* hit_v, int _id){

  pos_x = _pos_x;
  pos_y = _pos_y;
  d     = _d;
  p_x = hit_x;
  p_u = hit_u;
  p_v = hit_v;
  id  = _id;

}

FiberHitXUV::~FiberHitXUV() {}


FiberHitAna* FiberHitXUV::GetHit0(){
  double zx = p_x->GetZOrg();
  double zu = p_u->GetZOrg();
  double zv = p_v->GetZOrg();
  FiberHitAna* p_tmp = p_x;
  if(zx<zu && zx<zv) p_tmp = p_x;
  if(zu<zv && zu<zx) p_tmp = p_u;
  if(zv<zx && zv<zu) p_tmp = p_v;
  return p_tmp;
};

FiberHitAna* FiberHitXUV::GetHit1(){
  double zx = p_x->GetZOrg();
  double zu = p_u->GetZOrg();
  double zv = p_v->GetZOrg();
  FiberHitAna* p_tmp = p_u;
  if( (zx<zu && zx>zv) || (zx>zu && zx<zv) ) p_tmp = p_x;
  if( (zu<zv && zu>zx) || (zu>zv && zu<zx) ) p_tmp = p_u;
  if( (zv<zx && zv>zu) || (zv>zx && zv<zu) ) p_tmp = p_v;
  return p_tmp;
};

FiberHitAna* FiberHitXUV::GetHit2(){
  double zx = p_x->GetZOrg();
  double zu = p_u->GetZOrg();
  double zv = p_v->GetZOrg();
  FiberHitAna* p_tmp = p_v;
  if(zx>zu && zx>zv) p_tmp = p_x;
  if(zu>zv && zu>zx) p_tmp = p_u;
  if(zv>zx && zv>zu) p_tmp = p_v;
  return p_tmp;
};

FiberHitAna* FiberHitXUV::GetHit(int _i){
  FiberHitAna* p_tmp = p_x;
  if     (_i == 0) p_tmp = GetHit0();
  else if(_i == 1) p_tmp = GetHit1();
  else if(_i == 2) p_tmp = GetHit2();
  else p_tmp = nullptr;
  return p_tmp;
};

void FiberHitXUV::Print(void){
  std::cout << "-- HitXUV -----"  << std::endl;
  std::cout << "pos_x    : " << pos_x    << std::endl;
  std::cout << "pox_y    : " << pos_y    << std::endl;
  std::cout << "d        : " << d        << std::endl;
  p_x->Print();
  p_u->Print();
  p_v->Print();
  std::cout << "--------------"   << std::endl;
}
