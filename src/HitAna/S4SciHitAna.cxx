#include "ConstantParameter.hh"
#include "FiberHitAna.hh"
#include "TString.h"
#include <iostream>
#include <math.h>
#include "S4SciHitAna.hh"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

S4SciHitAna::S4SciHitAna(S4TQ *a, ParaManager *par){
  // tdc_sc31 = a.tdc_sc31; /// copy class members here
  // tdc_sc41 = a.tdc_sc41;
  // tdc_sc42 = a.tdc_sc42;
  // tdc_sc43 = a.tdc_sc43;
  // qdc_sc31 = a.qdc_sc31; /// copy class members here
  // qdc_sc41 = a.qdc_sc41;
  // qdc_sc42 = a.qdc_sc42;
  // qdc_sc42_lowgain = a.qdc_sc42_lowgain;
  // qdc_sc43 = a.qdc_sc43;
  // tdc_trig = a.tdc_trig;
  // tdc_s2tref = a.tdc_s2tref;

  t_sc31 = (par->mtdc_ch2ns)*0.5*(a->tdc_sc31[0]+a->tdc_sc31[1]);
  t_sc41 = (par->mtdc_ch2ns)*(0.25*((float)(a->tdc_sc41[0][0]+a->tdc_sc41[0][1]))
      +0.125*((float)(a->tdc_sc41[1][0]+a->tdc_sc41[1][1]+a->tdc_sc41[2][0]+a->tdc_sc41[2][1])));
  t_sc42 = (par->mtdc_ch2ns) * 0.5 * (a->tdc_sc42[0]+a->tdc_sc42[1]);
  t_sc43 = (par->mtdc_ch2ns) * 0.125 *((a->tdc_sc43[0][0]+a->tdc_sc43[0][1]) + (a->tdc_sc41[1][0]+a->tdc_sc41[1][1])
      + (a->tdc_sc43[2][0]+a->tdc_sc43[2][1]) + (a->tdc_sc41[3][0]+a->tdc_sc41[3][1]));

  de_sc31 = sqrt(a->qdc_sc31[0]*a->qdc_sc31[1]);
  de_sc41 = sqrt(a->qdc_sc41[0][0]*a->qdc_sc41[0][1]);
  //de_sc42_low  = sqrt(a->qdc_sc42_lowgain[0]*a->qdc_sc42_lowgain[1]);
  //de_sc42_high = sqrt(a->qdc_sc42[0]*a->qdc_sc42[1]);  // pedestal substraction here, or not needed?
  de_sc42_high  = sqrt(a->qdc_sc42_lowgain[0]*a->qdc_sc42_lowgain[1]);
  de_sc42_low   = sqrt(a->qdc_sc42[0]*a->qdc_sc42[1]);  // pedestal substraction here, or not needed?

  tof_sc3141 = t_sc41-t_sc31 + par->offset_tof_sc3141; ///
  tof_sc4143 = t_sc43-t_sc41 + par->offset_tof_sc4143; /// path length correction can be done in OpticsMom

  ResIdentify(); // get PID of residual w/ de_sc42 and tof_sc3141 cut
};

S4SciHitAna::~S4SciHitAna(){};

void S4SciHitAna::ResIdentify(){
  // TParticlePDG *partRoot = TDatabasePDG::Instance()->GetParticle(pdg[ipdg]);
  if(tof_sc3141>60 && tof_sc3141<65) {
    if(de_sc42_high>1500 && de_sc42_high<3500){  //He3 (H3L)
      PID_residual = 1000020030;
      mass_residual = 2808.391; // in MeV/(c^2)
    }
  }
  if(tof_sc3141>50 && tof_sc3141<60) {
    if(de_sc42_high>300 && de_sc42_high<1200){   // d (nnL)
      PID_residual = 1000010020;
      mass_residual = 1875.613;
    }
    else if(de_sc42_high>1500 && de_sc42_high<3500){ // He4 (H4L)
      PID_residual = 1000020040;
      mass_residual = 3727.379;
    }
  }
  if(tof_sc3141>50 && tof_sc3141<60) {
    if(de_sc42_high>400 && de_sc42_high<1500){  // p (Lambda)
      PID_residual = 2212;
      mass_residual = 938.27203; // in MeV/(c^2)
    }
  }

  /// Li6, C12, B9 w/ C-12 beam etc later..
};

void S4SciHitAna::Print(){
  std::cout << "--  S4SciHitAna  --"<< std::endl;
  std::cout << "tof_sc3141 : " << tof_sc3141 << std::endl;
  std::cout << "de_sc42_high : " << de_sc42_high << std::endl;
  std::cout << "PID : " << PID_residual << std::endl;
  std::cout << "-------------------"<< std::endl;

}
