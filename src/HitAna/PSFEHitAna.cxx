#include "ConstantParameter.hh"
#include "PSFEHitAna.hh"
#include "ParaManager.hh"
#include "TMath.h"
#include <iostream>
#include <math.h>

PSFEHitAna::PSFEHitAna(int seg, int t, int q, ParaManager *par, double t_t0){
  _seg     = seg;
  _t       = t;
  _q       = q;
  SetZRPhi(par);
  _time = t * par->psfe_ch2ns - t_t0 + par->wasa_tof_offset;
}

PSFEHitAna::~PSFEHitAna() {}


void PSFEHitAna::Print(void){
  std::cout << "--------------"  << std::endl;
  std::cout << "seg : " << _seg << std::endl;
  std::cout << "--------------"  << std::endl;
}

void PSFEHitAna::SetZRPhi(ParaManager *par){

  double phi_buf[44] = {
    4.09  , 12.27 , 20.45 , 28.63 , 36.81 , 44.99 , 53.17 , 61.35 , 69.53 , 77.71 ,
    85.89 , 94.07 , 102.25, 110.43, 118.61, 126.79, 134.97, 143.15, 151.33, 159.51,
    167.69, 175.87, 184.05, 192.23, 200.41, 208.59, 216.77, 224.95, 233.13, 241.31,
    249.49, 257.67, 265.85, 274.03, 282.21, 290.39, 298.57, 306.75, 314.93, 323.11,
    331.29, 339.47, 347.65, 355.83
  };

  _phi = phi_buf[_seg] * TMath::Pi() / 180.;
  if( _phi > TMath::Pi() ) _phi -= 2*TMath::Pi();
  _z = par->psfe_pos_z;
  _rmin = par->psfe_rmin;
  _rmax = par->psfe_rmax;

}
