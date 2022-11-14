#include "ConstantParameter.hh"
#include "PSBHitAna.hh"
#include "ParaManager.hh"
#include "TMath.h"
#include <iostream>
#include <math.h>

PSBHitAna::PSBHitAna(int seg, int t_u, int t_d, int q_u, int q_d, ParaManager *par){
  _seg     = seg;
  _t_u     = t_u;
  _t_d     = t_d;
  _q_u     = q_u;
  _q_d     = q_d;
  SetRPhi();
  _time_u = t_u * par->psb_ch2ns[seg][0] + par->psb_off_time[seg][0];
  _time_d = t_d * par->psb_ch2ns[seg][1] + par->psb_off_time[seg][1];
  _time   = (_time_u + _time_d)/2.;
  _z      = (_time_u - _time_d) * par->psb_zpar[seg];
  //_time = (t_u * par->psb_ch2ns + t_d * par->psb_ch2ns)/2. + par->psb_off_time[seg];
}

PSBHitAna::~PSBHitAna() {}


void PSBHitAna::Print(void){
  std::cout << "-- PSB ------------"  << std::endl;
  std::cout << "seg : " << _seg << std::endl;
  std::cout << "t_u : " << _t_u << std::endl;
  std::cout << "t_d : " << _t_d << std::endl;
  std::cout << "time_u : " << _time_u << std::endl;
  std::cout << "time_d : " << _time_d << std::endl;
  std::cout << "time   : " << _time   << std::endl;
  std::cout << "-------------------"  << std::endl;
}

void PSBHitAna::SetRPhi(void){

  if(_seg<23){
    if (0 == (_seg % 2)) { // inner PSB
      _r   = 221.25;
      _phi = TMath::Pi() * (9.15 + 14.7 * ((double)(_seg / 2))) / 180.0;
    } else { // outer PSB
      _r   = 231.75;
      _phi = TMath::Pi() * (16.5 + 14.7 * ((double)((_seg - 1) / 2))) / 180.0;
    }
  }
  else{
    if (0 == ((_seg-23) % 2)) { // inner PSB
      _r   =  221.25;
      _phi = TMath::Pi() * (189.15 + 14.7 * ((float)((_seg-23) / 2))) / 180.0;
    } else { // uter PSB
      _r   = 231.75;
      _phi = TMath::Pi() * (196.5 + 14.7 * ((float)(((_seg-23) - 1) / 2))) / 180.0;
    }
  }

  _phi += TMath::Pi()/2.;
  if( _phi > TMath::Pi() ) _phi -= 2*TMath::Pi();

}
