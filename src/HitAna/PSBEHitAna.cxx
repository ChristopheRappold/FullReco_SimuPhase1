#include "ConstantParameter.hh"
#include "PSBEHitAna.hh"
#include "ParaManager.hh"
#include "TMath.h"
#include <iostream>
#include <math.h>

PSBEHitAna::PSBEHitAna(int seg, int t, int q, ParaManager *par){
  _seg     = seg;
  _t       = t;
  _q       = q;
  SetZRPhi(par);
}

PSBEHitAna::~PSBEHitAna() {}


void PSBEHitAna::Print(void){
  std::cout << "--------------"  << std::endl;
  std::cout << "seg : " << _seg << std::endl;
  std::cout << "--------------"  << std::endl;
}

void PSBEHitAna::SetZRPhi(ParaManager *par){

  double phi_buf[38] = {
    0,      8.18,   16.36,  24.54,  32.72, 40.9,   49.08,  57.26,  65.44,  73.62,
    106.38, 114.56, 122.74, 130.92, 139.1, 147.28, 155.46, 163.64, 171.82, 180,
    188.18, 196.36, 204.54, 212.72, 220.9, 229.08, 237.26, 245.44, 253.62, 286.38,
    294.56, 302.74, 310.92, 319.1, 327.28, 335.46, 343.64, 351.82
  };

  _phi = phi_buf[_seg] * TMath::Pi() / 180.;
  if( _phi > TMath::Pi() ) _phi -= 2*TMath::Pi();
  _z = par->psbe_pos_z;
  _rmin = par->psbe_rmin;
  _rmax = par->psbe_rmax;

}
