#include "ConstantParameter.hh"
#include "FiberHitAna.hh"
#include "FiberHitXUV.hh"
#include "TString.h"
#include <iostream>
#include <math.h>

FiberHitXUV::FiberHitXUV(double _pos_x, double _pos_y, double _d, FiberHitAna* hit_x, FiberHitAna* hit_u, FiberHitAna* hit_v){

  pos_x = _pos_x;
  pos_y = _pos_y;
  d     = _d;
  p_x = hit_x;
  p_u = hit_u;
  p_v = hit_v;

}

FiberHitXUV::~FiberHitXUV() {}


FiberHitAna* FiberHitXUV::GetHit0(){
  double zx = p_x->GetZ();
  double zu = p_u->GetZ();
  double zv = p_v->GetZ();
  FiberHitAna* p_tmp = 0;
  if(zx<zu && zx<zv) p_tmp = p_x;
  if(zu<zv && zu<zx) p_tmp = p_u;
  if(zv<zx && zv<zu) p_tmp = p_v;
  return p_tmp;
};

FiberHitAna* FiberHitXUV::GetHit1(){
  double zx = p_x->GetZ();
  double zu = p_u->GetZ();
  double zv = p_v->GetZ();
  FiberHitAna* p_tmp = 0;
  if( (zx<zu && zx>zv) || (zx>zu && zx<zv) ) p_tmp = p_x;
  if( (zu<zv && zu>zx) || (zu>zv && zu<zx) ) p_tmp = p_u;
  if( (zv<zx && zv>zu) || (zv>zx && zv<zu) ) p_tmp = p_v;
  return p_tmp;
};

FiberHitAna* FiberHitXUV::GetHit2(){
  double zx = p_x->GetZ();
  double zu = p_u->GetZ();
  double zv = p_v->GetZ();
  FiberHitAna* p_tmp = 0;
  if(zx>zu && zx>zv) p_tmp = p_x;
  if(zu>zv && zu>zx) p_tmp = p_u;
  if(zv>zx && zv>zu) p_tmp = p_v;
  return p_tmp;
};

//void FiberHitXUV::Print(void){
//}
