#ifndef FIBER_HIT_ANA_HH
#define FIBER_HIT_ANA_HH
#include "../EventWASAUnpack/WASAUnpackBranch.hh"
#include "ParaManager.hh"
#include "limits.h"

class FiberHitAna
{
  public:
    FiberHitAna(FiberHit a, ParaManager *par, int tref, double t_t0);
    FiberHitAna(FiberHitAna *a);
    ~FiberHitAna();
    int    GetSfp(){   return i_sfp;};
    int    GetCtdc(){  return i_ctdc;};
    int    GetCh(){    return i_ch;};
    int    GetDet(){   return i_detector;};
    int    GetLay(){   return i_layer;};
    int    GetFib(){   return i_fiber;};
    int    GetClsize(){return i_clsize;};
    int    GetTL(){    return t_leading;};
    int    GetTT(){    return t_trailing;};
    double GetPos(){   return i_pos;};
    double GetAng(){   return i_ang;};
    double GetTOT(){   return t_tot;};
    double GetTime(){  return t_time;};
    double GetClFib(){ return i_cl_fiber;};
    double GetZ(){     return i_z;};
    double GetRes(){   return i_res;};
    int    GetDid(){   return i_did;};
    bool   IsValid(){  return i_valid;};
    void   Print();
    void   Add(FiberHitAna* hit);

  private:
    int i_sfp      = -1;
    int i_ctdc     = -1;
    int i_ch       = -1;
    int i_detector = -1;
    int i_layer    = -1;
    int i_fiber    = -1;
    int i_clsize   = 1;
    int t_leading  = INT_MIN;
    int t_trailing = INT_MIN;
    double i_pos = -9999.;
    double i_ang = -9999.;
    double t_tot = -9999.;
    double t_time = -9999.;
    double i_cl_fiber = -9999.;
    double i_z   = -9999.;
    double i_res = -9999.;
    int  i_did = -1;
    bool i_valid = false;
    double GetAngFiber(ParaManager *par);
    double GetPosFiber(ParaManager *par);
    double GetZFiber(  ParaManager *par);
    bool GetValidFiber(ParaManager *par);
};


#endif
