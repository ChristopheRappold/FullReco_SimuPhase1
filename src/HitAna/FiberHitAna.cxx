#include "ConstantParameter.hh"
#include "FiberHitAna.hh"
#include "TString.h"
#include <iostream>
#include <math.h>

FiberHitAna::FiberHitAna(FiberHit a, ParaManager *par, int tref, double t_t0){
  i_sfp      = a.i_sfp;
  i_ctdc     = a.i_ctdc;
  i_ch       = a.i_ch;
  i_detector = a.i_detector;
  i_layer    = a.i_layer;
  i_fiber    = a.i_fiber;
  i_clsize   = 1;
  t_leading  = a.t_leading  - tref;
  t_trailing = a.t_trailing - tref;
  i_ang = GetAngFiber(par);
  i_pos = GetPosFiber(par);
  t_tot = t_trailing - t_leading;
  t_time = (t_leading + par->fiber_time_offset[i_detector][i_layer][i_fiber]) * par->fiber_ch2ns - t_t0;
  //std::cout << t_leading << "\n";

  //std::cout << par->fiber_time_offset[i_detector][i_layer][i_fiber] << "\n";
  i_cl_fiber = a.i_fiber;
  i_z        = GetZFiber(par);
  //i_res      = 0.15;
  i_res      = par->fiber_res;
  i_did   = i_detector * 10 + i_layer;
  i_valid = GetValidFiber(par);
}

FiberHitAna::FiberHitAna(FiberHitAna *a){
  i_sfp      = a->GetSfp();
  i_ctdc     = a->GetCtdc();
  i_ch       = a->GetCh();
  i_detector = a->GetDet();
  i_layer    = a->GetLay();
  i_fiber    = a->GetFib();
  i_clsize   = a->GetClsize();
  t_leading  = a->GetTL();
  t_trailing = a->GetTT();
  i_ang      = a->GetAng();
  i_pos      = a->GetPos();
  t_tot      = a->GetTOT();
  t_time     = a->GetTime();
  i_cl_fiber = a->GetClFib();
  i_z        = a->GetZ();
  i_res      = a->GetRes();
  i_did      = a->GetDid();
  i_valid    = a->IsValid();
}

FiberHitAna::~FiberHitAna() {}


void FiberHitAna::Print(void){
  std::cout << "- Fiber ------"  << std::endl;
  std::cout << "det : " << i_detector << std::endl;
  std::cout << "lay : " << i_layer << std::endl;
  std::cout << "fib : " << i_fiber << std::endl;
  std::cout << "clsize : " << i_clsize << std::endl;
  std::cout << Form("fib_cl : %.1f", i_cl_fiber)  << std::endl;
  std::cout << Form("pos    : %.1f", i_pos)  << std::endl;
  std::cout << Form("z      : %.1f", i_z)  << std::endl;
  std::cout << "--------------"  << std::endl;
}

double FiberHitAna::GetPosFiber(ParaManager *par){

  double pos = -9999.;
  double step = 1.10;

  switch(i_detector){

    // UFT1
    case 0:
      if(i_fiber%2==0) pos = -0.1375 - 0.275*2 + (i_fiber/2 - 128/2 +1) *1.1;
      else             pos = -0.1375 - 0.275*1 + (i_fiber/2 - 128/2 +1) *1.1;
      switch(i_layer){
        case 0: pos =  pos + par->fiber_uft1_pos_x*cos(i_ang) + par->fiber_uft1_pos_y*sin(i_ang) + par->fiber_uft1_off_x; break;
        case 1: pos =  pos + par->fiber_uft1_pos_x*cos(i_ang) + par->fiber_uft1_pos_y*sin(i_ang) + par->fiber_uft1_off_u; break;
        case 2: pos =  pos + par->fiber_uft1_pos_x*cos(i_ang) + par->fiber_uft1_pos_y*sin(i_ang) + par->fiber_uft1_off_v; break;
        default: break;
      }
      pos = pos + par->fiber_offset[i_detector][i_layer][i_fiber];
      break;

    // UFT2
    case 1:
      if(i_fiber%2==0) pos = -0.1375 - 0.275*2 + (i_fiber/2 - 128/2 +1) *1.1;
      else             pos = -0.1375 - 0.275*1 + (i_fiber/2 - 128/2 +1) *1.1;
      switch(i_layer){
        case 0: pos =  pos + par->fiber_uft2_pos_x*cos(i_ang) + par->fiber_uft2_pos_y*sin(i_ang) + par->fiber_uft2_off_x; break;
        case 1: pos =  pos + par->fiber_uft2_pos_x*cos(i_ang) + par->fiber_uft2_pos_y*sin(i_ang) + par->fiber_uft2_off_u; break;
        case 2: pos =  pos + par->fiber_uft2_pos_x*cos(i_ang) + par->fiber_uft2_pos_y*sin(i_ang) + par->fiber_uft2_off_v; break;
        default: break;
      }
      pos = pos + par->fiber_offset[i_detector][i_layer][i_fiber];
      break;

    // UFT3
    case 2:
      if(i_fiber%2==0) pos = 0.1375 + 0.275*2 - (i_fiber/2 - 192/2 +1) *1.1;
      else             pos = 0.1375 + 0.275*1 - (i_fiber/2 - 192/2 +1) *1.1;
      switch(i_layer){
        case 0: pos =  pos + par->fiber_uft3_pos_x*cos(i_ang) + par->fiber_uft3_pos_y*sin(i_ang) + par->fiber_uft3_off_x; break;
        case 1: pos =  pos + par->fiber_uft3_pos_x*cos(i_ang) + par->fiber_uft3_pos_y*sin(i_ang) + par->fiber_uft3_off_u; break;
        case 2: pos =  pos + par->fiber_uft3_pos_x*cos(i_ang) + par->fiber_uft3_pos_y*sin(i_ang) + par->fiber_uft3_off_v; break;
        default: break;
      }
      pos = pos + par->fiber_offset[i_detector][i_layer][i_fiber];
      break;

    // MFT1,2
    case 3:
    case 4:
      switch((i_detector-3)*6 +  i_layer*2 + (i_fiber>=128) ){
        case  0: step =  par->fiber_mft1_step_x1; break;
        case  1: step =  par->fiber_mft1_step_x2; break;
        case  2: step =  par->fiber_mft1_step_u1; break;
        case  3: step =  par->fiber_mft1_step_u2; break;
        case  4: step =  par->fiber_mft1_step_v1; break;
        case  5: step =  par->fiber_mft1_step_v2; break;
        case  6: step =  par->fiber_mft2_step_x1; break;
        case  7: step =  par->fiber_mft2_step_x2; break;
        case  8: step =  par->fiber_mft2_step_v1; break;
        case  9: step =  par->fiber_mft2_step_v2; break;
        case 10: step =  par->fiber_mft2_step_u1; break;
        case 11: step =  par->fiber_mft2_step_u2; break;
        default: break;
      }
      if(i_fiber< 128){
        if(i_fiber%2==0) pos = -(12.3 + 0.25 + 1.1/4*2) - (63 - i_fiber/2) *step;
        else             pos = -(12.3 + 0.25 + 1.1/4*1) - (63 - i_fiber/2) *step;
      }
      else{
        if(i_fiber%2==0) pos =  (12.3 + 0.25 + 1.1/4*1) + (i_fiber/2 - 64) *step;
        else             pos =  (12.3 + 0.25 + 1.1/4*2) + (i_fiber/2 - 64) *step;
      }
      switch((i_detector-3)*6 +  i_layer*2 + (pos>0) ){
        case  0: pos = pos + par->fiber_mft1_pos_x*cos(i_ang) + par->fiber_mft1_pos_y*sin(i_ang) + par->fiber_mft1_off_x1; break;
        case  1: pos = pos + par->fiber_mft1_pos_x*cos(i_ang) + par->fiber_mft1_pos_y*sin(i_ang) + par->fiber_mft1_off_x2; break;
        case  2: pos = pos + par->fiber_mft1_pos_x*cos(i_ang) + par->fiber_mft1_pos_y*sin(i_ang) + par->fiber_mft1_off_u1; break;
        case  3: pos = pos + par->fiber_mft1_pos_x*cos(i_ang) + par->fiber_mft1_pos_y*sin(i_ang) + par->fiber_mft1_off_u2; break;
        case  4: pos = pos + par->fiber_mft1_pos_x*cos(i_ang) + par->fiber_mft1_pos_y*sin(i_ang) + par->fiber_mft1_off_v1; break;
        case  5: pos = pos + par->fiber_mft1_pos_x*cos(i_ang) + par->fiber_mft1_pos_y*sin(i_ang) + par->fiber_mft1_off_v2; break;
        case  6: pos = pos + par->fiber_mft2_pos_x*cos(i_ang) + par->fiber_mft2_pos_y*sin(i_ang) + par->fiber_mft2_off_x1; break;
        case  7: pos = pos + par->fiber_mft2_pos_x*cos(i_ang) + par->fiber_mft2_pos_y*sin(i_ang) + par->fiber_mft2_off_x2; break;
        case  8: pos = pos + par->fiber_mft2_pos_x*cos(i_ang) + par->fiber_mft2_pos_y*sin(i_ang) + par->fiber_mft2_off_v1; break;
        case  9: pos = pos + par->fiber_mft2_pos_x*cos(i_ang) + par->fiber_mft2_pos_y*sin(i_ang) + par->fiber_mft2_off_v2; break;
        case 10: pos = pos + par->fiber_mft2_pos_x*cos(i_ang) + par->fiber_mft2_pos_y*sin(i_ang) + par->fiber_mft2_off_u1; break;
        case 11: pos = pos + par->fiber_mft2_pos_x*cos(i_ang) + par->fiber_mft2_pos_y*sin(i_ang) + par->fiber_mft2_off_u2; break;
        default: break;
      }
      pos = pos + par->fiber_offset[i_detector][i_layer][i_fiber];
      break;

    // DFT1
    case 5:
      if(i_fiber%2==0) pos =  0.1375 + 0.275*2 - (i_fiber/2 - 128/2 +1) *1.1;
      else             pos =  0.1375 + 0.275*1 - (i_fiber/2 - 128/2 +1) *1.1;
      switch(i_layer){
        case 0: pos =  pos + par->fiber_dft1_pos_x*cos(i_ang) + par->fiber_dft1_pos_y*sin(i_ang) + par->fiber_dft1_off_v; break;
        case 1: pos =  pos + par->fiber_dft1_pos_x*cos(i_ang) + par->fiber_dft1_pos_y*sin(i_ang) + par->fiber_dft1_off_u; break;
        case 2: pos =  pos + par->fiber_dft1_pos_x*cos(i_ang) + par->fiber_dft1_pos_y*sin(i_ang) + par->fiber_dft1_off_x; break;
        default: break;
      }
      pos = pos + par->fiber_offset[i_detector][i_layer][i_fiber];
      break;

    // DFT2
    case 6:
      if(i_fiber%2==0) pos =  0.1375 + 0.275*2 - (i_fiber/2 - 128/2 +1) *1.1;
      else             pos =  0.1375 + 0.275*1 - (i_fiber/2 - 128/2 +1) *1.1;
      switch(i_layer){
        case 0: pos =  pos + par->fiber_dft2_pos_x*cos(i_ang) + par->fiber_dft2_pos_y*sin(i_ang) + par->fiber_dft2_off_x; break;
        case 1: pos =  pos + par->fiber_dft2_pos_x*cos(i_ang) + par->fiber_dft2_pos_y*sin(i_ang) + par->fiber_dft2_off_u; break;
        case 2: pos =  pos + par->fiber_dft2_pos_x*cos(i_ang) + par->fiber_dft2_pos_y*sin(i_ang) + par->fiber_dft2_off_v; break;
        default: break;
      }
      pos = pos + par->fiber_offset[i_detector][i_layer][i_fiber];
      break;

  }

  return pos;

}

double FiberHitAna::GetAngFiber(ParaManager *par){

  double ang[7][3] = {
    {  0,  30, -30}, /* UFT1 */
    {  0,  30, -30}, /* UFT2 */
    {  0,  30, -30}, /* UFT3 */
    {  0, -60,  60}, /* MFT1 */
    {  0,  60, -60}, /* MFT2 */
    { 30, -30,   0}, /* DFT1 */
    {  0,  30, -30}  /* DFT2 */
  };

  double ang_offset = 0;
  if(i_detector==3 || i_detector==4 ){
    switch((i_detector-3)*6 +  i_layer*2 + (i_fiber>=128) ){
      case  0: ang_offset =  par->fiber_mft1_off_ang_x1; break;
      case  1: ang_offset =  par->fiber_mft1_off_ang_x2; break;
      case  2: ang_offset =  par->fiber_mft1_off_ang_u1; break;
      case  3: ang_offset =  par->fiber_mft1_off_ang_u2; break;
      case  4: ang_offset =  par->fiber_mft1_off_ang_v1; break;
      case  5: ang_offset =  par->fiber_mft1_off_ang_v2; break;
      case  6: ang_offset =  par->fiber_mft2_off_ang_x1; break;
      case  7: ang_offset =  par->fiber_mft2_off_ang_x2; break;
      case  8: ang_offset =  par->fiber_mft2_off_ang_v1; break;
      case  9: ang_offset =  par->fiber_mft2_off_ang_v2; break;
      case 10: ang_offset =  par->fiber_mft2_off_ang_u1; break;
      case 11: ang_offset =  par->fiber_mft2_off_ang_u2; break;
      default: break;
    }
  }

  return (ang[i_detector][i_layer] + ang_offset) * Deg2Rad;

}

double FiberHitAna::GetZFiber(ParaManager *par){

  double pos_z = -9999.;
  switch(i_detector){
    case 0: pos_z = par->fiber_uft1_pos_z + 4*(i_layer-1); break;
    case 1: pos_z = par->fiber_uft2_pos_z + 4*(i_layer-1); break;
    case 2: pos_z = par->fiber_uft3_pos_z + 4*(i_layer-1); break;
    case 3: pos_z = par->fiber_mft1_pos_z + 4*(i_layer-1); break;
    case 4: pos_z = par->fiber_mft2_pos_z + 4*(i_layer-1); break;
    case 5: pos_z = par->fiber_dft1_pos_z + 4*(i_layer-1); break;
    case 6: pos_z = par->fiber_dft2_pos_z + 4*(i_layer-1); break;
    default: break;
  }

  return pos_z;

}

bool FiberHitAna::GetValidFiber(ParaManager *par){
  bool judge = true;

  switch(i_sfp){
    case 0: if( par->fiber_mft_cut_min   > t_leading || par->fiber_mft_cut_max   < t_leading) judge = false; break;
    case 1: if( par->fiber_uft3_cut_min  > t_leading || par->fiber_uft3_cut_max  < t_leading) judge = false; break;
    case 2: if( par->fiber_uft12_cut_min > t_leading || par->fiber_uft12_cut_max < t_leading) judge = false; break;
    case 3: if( par->fiber_dft12_cut_min > t_leading || par->fiber_dft12_cut_max < t_leading) judge = false; break;
    default: break;
  }

  switch(i_sfp){
    case 0: if( par->fiber_tot_mft_cut_min   > t_tot || par->fiber_tot_mft_cut_max   < t_tot) judge = false; break;
    case 1: if( par->fiber_tot_uft3_cut_min  > t_tot || par->fiber_tot_uft3_cut_max  < t_tot) judge = false; break;
    case 2: if( par->fiber_tot_uft12_cut_min > t_tot || par->fiber_tot_uft12_cut_max < t_tot) judge = false; break;
    case 3: if( par->fiber_tot_dft12_cut_min > t_tot || par->fiber_tot_dft12_cut_max < t_tot) judge = false; break;
    default: break;
  }

  switch(i_sfp){
    case 0: if( par->fiber_time_mft_cut_min   > t_time || par->fiber_time_mft_cut_max   < t_time) judge = false; break;
    case 1: if( par->fiber_time_uft3_cut_min  > t_time || par->fiber_time_uft3_cut_max  < t_time) judge = false; break;
    case 2: if( par->fiber_time_uft12_cut_min > t_time || par->fiber_time_uft12_cut_max < t_time) judge = false; break;
    case 3: if( par->fiber_time_dft12_cut_min > t_time || par->fiber_time_dft12_cut_max < t_time) judge = false; break;
    default: break;
  }

  return judge;
}

void FiberHitAna::Add(FiberHitAna *hit){

  if(t_tot < hit->GetTOT()){
    i_ch    = hit->GetCh();
    i_fiber = hit->GetFib();
    t_leading  = hit->GetTL();
    t_trailing = hit->GetTT();
    t_tot    = hit->GetTOT();
  }

  i_pos      = ( i_pos      * i_clsize + hit->GetPos()   * hit->GetClsize() ) / ( i_clsize + hit->GetClsize() );
  i_cl_fiber = ( i_cl_fiber * i_clsize + hit->GetClFib() * hit->GetClsize() ) / ( i_clsize + hit->GetClsize() );
  i_clsize += hit->GetClsize();

}

