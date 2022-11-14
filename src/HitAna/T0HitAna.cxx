#include "ConstantParameter.hh"
#include "T0HitAna.hh"
#include "ParaManager.hh"
#include "TMath.h"
#include <iostream>
#include <math.h>

T0HitAna::T0HitAna(int seg, int t_u, int t_d, int q_u, int q_d, ParaManager *par){
  _seg     = seg;
  _t_u     = t_u;
  _t_d     = t_d;
  _q_u     = q_u;
  _q_d     = q_d;
  _z       = par->t0_pos_z;
  _time_u = t_u * par->t0_ch2ns[seg][0] + par->t0_off_time[seg][0];
  _time_d = t_d * par->t0_ch2ns[seg][1] + par->t0_off_time[seg][1];
  _time   = (_time_u + _time_d)/2.;
}

T0HitAna::~T0HitAna() {}


void T0HitAna::Print(void){
  std::cout << "--------------"  << std::endl;
  std::cout << "seg : " << _seg << std::endl;
  std::cout << "--------------"  << std::endl;
}
