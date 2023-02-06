#include "ConstantParameter.hh"
#include "TString.h"
#include <iostream>
#include <math.h>
#include "OpticsMom.hh"

OpticsMom::OpticsMom(FiberTrackAna *fiber, MWDCTracking  *mwdc, S4SciHitAna *s4scin, ParaManager *par){
  /// --- Optics Parameters ---- 17 Tm ---- ///
  par_optics["xx"] = par->optics_par["xx"];
  par_optics["xa"] = par->optics_par["xa"];
  par_optics["ax"] = par->optics_par["ax"];
  par_optics["aa"] = par->optics_par["aa"];
  par_optics["xd"] = par->optics_par["xd"];
  par_optics["ad"] = par->optics_par["ad"];
  par_optics["by"] = par->optics_par["by"];
  par_optics["bb"] = par->optics_par["bb"];
  par_optics["bd"] = par->optics_par["bd"];
  par_optics["xxd"] = par->optics_par["xxd"];
  par_optics["xad"] = par->optics_par["xad"];
  par_optics["axd"] = par->optics_par["axd"];
  par_optics["aad"] = par->optics_par["aad"];
  par_optics["byd"] = par->optics_par["byd"];
  par_optics["bbd"] = par->optics_par["bbd"];
  par_optics["xxxd"] = par->optics_par["xxxd"];
  par_optics["xaad"] = par->optics_par["xaad"];
  par_optics["axxd"] = par->optics_par["axxd"];
  par_optics["aaad"] = par->optics_par["aaad"];
  par_optics["bbbd"] = par->optics_par["bbbd"];
  par_optics["yy"] = par->optics_par["yy"];
  par_optics["yb"] = par->optics_par["yb"];
  par_optics["yd"] = par->optics_par["yd"];
  //// ---------
  //float X4_offset =  0.13;
  //float A4_offset =  -0.47 + 0.15;
  //float Y4_offset = -2.89;
  //float B4_offset = -2.01;
  float X4_offset = 0.38;
  float A4_offset = -0.47;
  float Y4_offset = -3.1;
  float B4_offset = -1.7;

  X4 = mwdc->GetX() + mwdc->GetA()*(par->dist_focS4 - par->dist_MWDC_zref) +  par->mwdc_shift_x_alignment;;
  A4 = mwdc->GetA()*1000. + par->mwdc_shift_a_alignment;
  Y4 = mwdc->GetY() + mwdc->GetB()*(par->dist_focS4 - par->dist_MWDC_zref) +  par->mwdc_shift_y_alignment;;
  B4 = mwdc->GetB()*1000. + par->mwdc_shift_b_alignment;
  chi2_s4 = mwdc->GetChi2();

  //// offset by limit S2 phase space, just to compensate the residual of reconstucted delta, A2, B2.
  X4 += X4_offset;
  A4 += A4_offset;
  Y4 += Y4_offset;
  B4 += B4_offset;

  X2 = fiber->GetX()+fiber->GetA()*par->optics_s2z; //par->dft12_z
  A2 = fiber->GetA()*1000.;
  Y2 = fiber->GetY()+fiber->GetB()*par->optics_s2z;
  B2 = fiber->GetB() *1000.;
  Z2 = par->optics_s2z;
  chi2_dft = fiber->GetChi2();
  tof_sc3141 = s4scin->GetTOF_sc3141();
  CalDelta();
  CorrectTOF();
  PID = s4scin->GetPID();
  mass = s4scin->GetMass() ; //He3 mass
  GetPref();
  CalMom();
  CalVec();
  track = fiber;
};
OpticsMom::~OpticsMom(){};

void OpticsMom::Print(){
  std::cout << "-- Optics ---------"  << std::endl;
  std::cout << "delta : " << delta << std::endl;
  std::cout << "-------------------"  << std::endl;
};

void OpticsMom::CalDelta(){
  float coeff_delta = 0.37; /// from DFT to S4, we use 0.4, from vacuum window to S4, we use 0.37
  float a_tmp = par_optics["xd"][1]
    // + par_optics["xxxd"][1] * pow(A4 - coeff_delta * X2, 2)
    // + par_optics["xaad"][1] * pow(A4, 2)
    + par_optics["xxd"][1] * (A4 - coeff_delta * X2) + par_optics["xad"][1] * A4;
  float b_tmp = par_optics["xd"][0]
    // + par_optics["xxxd"][0] * pow(A4 - coeff_delta * X2, 2)
    // + par_optics["xaad"][0] * pow(A4, 2)
    + par_optics["xxd"][0] * (A4 - coeff_delta * X2) + par_optics["xad"][0] * A4;
  float c_tmp = par_optics["xx"][0] * (A4 - coeff_delta * X2) + par_optics["xx"][1] * pow(A4 - coeff_delta * X2, 2)
    + par_optics["xa"][0] * A4 + par_optics["xa"][1] * pow(A4, 2) - X4;
  float delta_tmp = (-1. / (2 * a_tmp)) * (b_tmp - sqrt(b_tmp * b_tmp - 4 * a_tmp * c_tmp));
  // to correct the energy loss difference in SC31 between He4 and d
  if(PID==1000020020) delta_tmp-=(0.04393+0.00072*delta_tmp)/100.;
  delta = delta_tmp;
};

void OpticsMom::GetPref(){
  if(     PID == 2212)       Pref = 2000.0; // Lambda(p), Brho(S2->S3)=6.6713 Tm, in MeV/c
  else if(PID == 1000020030) Pref = 7501.2; // H3L(He3), Brho(S2->S3)=12.51072 Tm, in MeV/c
  else if(PID == 1000020040) Pref = 9971.0; // H4L(He4), Brho(S2->S3)=16.6298 Tm, in MeV/c
  else if(PID == 1000020020) Pref = 4985.5; // nnL(d), Brho(S2->S3)=16.6298 Tm, in MeV/c
  else {
    std::cout<<" NO corresponding Pref found !!!!"<<std::endl;
    std::cout << "PID : " << PID << std::endl;
    Pref = 0.;
  }
}

void OpticsMom::CalMom(){
  if(Pref!=0.)
    momentum = (1. + delta) * Pref /1000.; // in GeV/c
  else{
    std::cout<<"Pref doesn't exsit, please check analys!!!"<<std::endl;
    momentum = 0.;
  }
}

void OpticsMom::CalVec(){
  double a = GetA2Rec()/1000.;
  double b = GetB2Rec()/1000.;
  double c = sqrt(1 + pow(a,2.) + pow(b,2.));
  double px = momentum * a/c;
  double py = momentum * b/c;
  double pz = momentum * 1/c;
  v_mom = TVector3(px, py, pz);
}

float OpticsMom::GetA2Rec(){
  float a_tmp = par_optics["aa"][1]
    + par_optics["aaad"][0] * delta + par_optics["aaad"][1] * pow(delta, 2)
    + par_optics["axxd"][0] * delta + par_optics["axxd"][1] * pow(delta, 2)
    + par_optics["ax"][1] ;
  float b_tmp = par_optics["aa"][0] + par_optics["ax"][0] - 2.0 * 0.45 * par_optics["ax"][1] * X2
    + par_optics["aad"][0] * delta + par_optics["aad"][1] * pow(delta, 2)
    - 2.0*0.45*X2*par_optics["axxd"][0]*delta
    - 2.0*0.45*X2*par_optics["axxd"][1]*pow(delta,2)
    + par_optics["axd"][0] * delta + par_optics["axd"][1] * pow(delta, 2);
  float c_tmp = par_optics["ax"][1] * pow(0.45 * X2, 2) - par_optics["ax"][0] * X2 * 0.45
    - 0.45*par_optics["axd"][0]*X2*delta-0.45*par_optics["axd"][1]*X2*pow(delta,2)
    + pow(0.45 * X2, 2) * par_optics["axxd"][0] * delta
    + pow(0.45 * X2, 2) * par_optics["axxd"][1] * pow(delta, 2)
    - A4 + par_optics["ad"][0] * delta + par_optics["ad"][1] * pow(delta, 2);
  float A2_recon = (-1.0 / (2.0 * a_tmp)) * (b_tmp - sqrt(b_tmp * b_tmp - 4.0 * a_tmp * c_tmp));
  return A2_recon;
}

float OpticsMom::GetB2Rec(){
  float coeff_b = 2.3; // DFT->S4: 2.0;  Vacuum window->S4: 2.3
  float B2_recon = 1./coeff_b* (Y2 - par_optics["by"][0]*Y4 - par_optics["bb"][0]*B4 - par_optics["bb"][1]*pow(B4,2)
      - par_optics["byd"][0]*Y4*delta - par_optics["byd"][1]*Y4*pow(delta,2)
      - par_optics["bbd"][0]*B4*delta - par_optics["bbd"][1]*B4*pow(delta,2)
      //  - par_optics["bbbd"][0]*delta*pow(B4,2) - par_optics["bbbd"][1]*pow(delta,2)*pow(B4,2)
      - par_optics["bd"][0]*delta - par_optics["bd"][1]*pow(delta,2));
  /// currently only 1 st order for B2 recon
  return B2_recon;
}

float OpticsMom::GetY2Rec(){
  float s4b_tmp = B4+1.3*Y4;
  float s4y_tmp = Y4;
  float Y2_recon = par_optics["yd"][0]*delta + par_optics["yd"][1]*pow(delta,2)
    + par_optics["yy"][0]*s4y_tmp + par_optics["yy"][1]*pow(s4y_tmp,2)
    + par_optics["yb"][0]*s4b_tmp + par_optics["yb"][1]*pow(s4b_tmp,2);
  return Y2_recon;
}

void OpticsMom::CorrectTOF(){ //// not yet completed
  double tof_corr_tmp = 0.;
  ////------
  tof_corr = tof_sc3141 - tof_corr_tmp;
}
