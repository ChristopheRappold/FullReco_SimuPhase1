#include "ConstantParameter.hh"
#include "MDCHitAna.hh"
#include "TString.h"
#include "TMath.h"
#include <iostream>
#include <math.h>

MDCHitAna::MDCHitAna(MDCHit a, ParaManager *par, int tref, double t_t0){
  i_ctdc     = a.i_ctdc;
  i_ch       = a.ch_ctdc;
  i_layer    = a.i_layer;
  i_wire     = a.i_wire;
  t_leading  = a.t_leading  - tref;
  t_trailing = a.t_trailing - tref;
  t_tot      = a.t_trailing - a.t_leading;
  i_valid = GetValidMDC(par);
  SetPhys(par);
  SetDriftTime(par, t_t0);
  SetDriftLength(par);
  //i_dl = i_wsize/2.;
  //i_dl = 0;
  //i_res = 0.2;
  //i_res = i_wsize*2/sqrt(12.);
  //i_res = 0.5;
  i_res = par->mdc_res;
  i_did = 70 + i_layer;
}

MDCHitAna::~MDCHitAna() {}


void MDCHitAna::Print(void){
  std::cout << "- MDC --------"  << std::endl;
  std::cout << "lay : " << i_layer << std::endl;
  std::cout << "wir : " << i_wire  << std::endl;
  std::cout << "r   : " << i_r     << std::endl;
  std::cout << "phi : " << i_phi   << std::endl;
  std::cout << "dphi f : " << i_dpf   << std::endl;
  std::cout << "dphi b : " << i_dpb   << std::endl;
  std::cout << "lf     : " << i_lf   << std::endl;
  std::cout << "lb     : " << i_lb   << std::endl;
  std::cout << "--------------"  << std::endl;
}

bool MDCHitAna::GetValidMDC(ParaManager *par){

  bool judge = true;
  switch(i_layer){
    case 0 : case 1 : case 2 : case 3 : case 4 :          if( par->mdc_cut_min > t_leading || par->mdc_cut_max1 < t_leading) judge = false; break;
    case 5 : case 6 : case 7 : case 8 : case 9 : case 10: if( par->mdc_cut_min > t_leading || par->mdc_cut_max2 < t_leading) judge = false; break;
    case 11: case 12: case 13: case 14: case 15: case 16: if( par->mdc_cut_min > t_leading || par->mdc_cut_max3 < t_leading) judge = false; break;
    default: break;
  }
  if( par->mdc_cut_tot > t_tot ) judge = false;

  return judge;
}


void MDCHitAna::SetPhys(ParaManager *par){

  const auto& Physinfo = par->mdc_PhysMDC[i_layer];

  i_z = par->mdc_pos_z;

  int num_wire = Physinfo.nb_wire;
  double rmax_wire = Physinfo.rmax;
  double size_wire = 0.5 * Physinfo.wire_size;//convert to radius
  double theta_begin_gr1 = 0.5*TMath::Pi() + hitana::Deg2Rad*Physinfo.Dpsi1; // or Dphi1??
  double theta_begin_gr2 = theta_begin_gr1 + TMath::Pi();
  double theta_step = hitana::Deg2Rad * Physinfo.dpsi0; // or dphi0??

  i_r = rmax_wire;

  double theta_tmp;
  if(i_wire<num_wire/2){ //Group 1 of this layer
    theta_tmp = theta_begin_gr1 + ((double)i_wire) * theta_step;
  }else{ // Group 2 of this layer
    theta_tmp = theta_begin_gr2 + ((double)(i_wire-(num_wire/2))) * theta_step;
  }
  if(theta_tmp>2.*TMath::Pi()) theta_tmp -= 2.*TMath::Pi();
  if(theta_tmp>TMath::Pi()) theta_tmp -= 2.*TMath::Pi();
  if(theta_tmp>0) theta_tmp =  TMath::Pi() - theta_tmp;
  else            theta_tmp = -TMath::Pi() - theta_tmp;

  i_phi = theta_tmp;

  i_dpf   = -Physinfo.Dpsi0f * hitana::Deg2Rad;
  i_dpb   = -Physinfo.Dpsi0b * hitana::Deg2Rad;
  i_lf    = -Physinfo.l_front;
  i_lb    = -Physinfo.l_back;
  i_zf    = -Physinfo.z_front;
  i_zb    = -Physinfo.z_back;
  i_wsize =  Physinfo.wire_size / 2.;

}

void MDCHitAna::SetDriftTime(ParaManager *par, double t_t0){

  double drift_time = t_leading * par->mdc_ch2ns -t_t0 - par->mdc_t0_off[i_layer] - par->mdc_t0_off_wir[i_layer][i_wire];

  i_dt =  drift_time;

}

void MDCHitAna::SetDriftLength(ParaManager *par){

  double par_dl[7];
  par_dl[0] = par->mdc_PhysMDC[i_layer].drift_par0;
  par_dl[1] = par->mdc_PhysMDC[i_layer].drift_par1;
  par_dl[2] = par->mdc_PhysMDC[i_layer].drift_par2;
  par_dl[3] = par->mdc_PhysMDC[i_layer].drift_par3;
  par_dl[4] = par->mdc_PhysMDC[i_layer].drift_par4;
  par_dl[5] = par->mdc_PhysMDC[i_layer].drift_par5;
  par_dl[6] = par->mdc_PhysMDC[i_layer].drift_par6;

  //std::cout << "par_dl : " << par_dl[1] << std::endl;
  //std::cout << "i_dt : " << i_dt << std::endl;

  double dl_buf = 0;
  for(int i=0; i<6; ++i){
    dl_buf += par_dl[i] * pow(i_dt - par_dl[6], i);
  }
  if(dl_buf<0)       dl_buf = 0;
  if(dl_buf>i_wsize) dl_buf = i_wsize;

  i_dl = dl_buf;

}
