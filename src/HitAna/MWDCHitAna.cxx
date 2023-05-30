#include "ConstantParameter.hh"
#include "MWDCHitAna.hh"
#include "TString.h"
#include "TMath.h"
#include <iostream>
#include <math.h>

MWDCHitAna::MWDCHitAna(MWDCHit a, ParaManager *par, int tref){
  i_ctdc = a.i_ctdc;
  i_ch = a.ch_ctdc;
  i_plane= a.i_plane;
  i_wire = a.i_wire;
  t_leading = a.t_leading - tref;
  t_trailing = a.t_trailing - tref;
  tot = t_trailing - t_leading;
  SetDriftTime(par);
  SetDriftLength(par);
  /// add paramters of each plane also here?
};
MWDCHitAna::~MWDCHitAna(){};

void MWDCHitAna::SetDriftTime(ParaManager *par){
  double drift_time = t_leading - par->mwdc_t0_off[i_plane]; 
  /// here we use t_leading w/o offset and ch2ns...
  i_dt = drift_time;
};

void MWDCHitAna::SetDriftLength(ParaManager *par){
  int npar = 8;
  double par_dl[16][npar]; // load from txt in ParaManager
  for(int i = 0; i < npar ; ++i)
    par_dl[i_plane][i] = par->mwdc_dtdx_par[i_plane][i];

  double dl_tmp = 0.;
  for(int i = 0; i < npar; ++i)
    dl_tmp += par_dl[i_plane][i]*pow(i_dt,i);

  if(dl_tmp<0) dl_tmp = 0.;
  if(dl_tmp>2.5) dl_tmp = 2.5; //max_dl
  i_dl = dl_tmp;
  i_dl = fabs(i_dl);
};
