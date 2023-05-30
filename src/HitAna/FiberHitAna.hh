#ifndef FIBER_HIT_ANA_HH
#define FIBER_HIT_ANA_HH
#include "../EventWASAUnpack/WASAUnpackBranch.hh"
#include "ParaManager.hh"
#include "limits.h"

class FiberHitAna
{
  public:
    FiberHitAna(FiberHit a, ParaManager *par, int tref, double t_t0 = 0);
    FiberHitAna(FiberHitAna *a);
    FiberHitAna(int detector, int layer, int fiber, double time, double energy, int simtrackid, ParaManager *par);
    FiberHitAna(double pos, double ang, double z, double res);
    ~FiberHitAna();
    int    GetSfp(){   return i_sfp;};
    int    GetCtdc(){  return i_ctdc;};
    int    GetCh(){    return i_ch;};
    int    GetDet(){   return i_detector;};
    int    GetLay(){   return i_layer;};
    int    GetFib(){   return i_fiber;};
    int    GetClsize(){return i_clsize;};
    int    GetUD(){    return i_ud;};
    int    GetTL(){    return t_leading;};
    int    GetTT(){    return t_trailing;};
    double GetPos(){   return i_pos;};
    double GetPosOrg(){return i_pos_org;};
    double GetPosL(){  return i_posl;};
    double GetdE(){    return i_dE;};
    double GetAng(){   return i_ang;};
    double GetTOT(){   return t_tot;};
    double GetTime(){  return t_time;};
    double GetClFib(){ return i_cl_fiber;};
    double GetZ(){     return i_z;};
    double GetZOrg(){  return i_z_org;};
    double GetRes(){   return i_res;};
    int    GetDid(){   return i_did;};
    bool   IsValid(){  return i_valid;};
    bool   IsUsed(){   return i_used;};
    int    GetSimTrackID(){return i_simtrackid;};
    void   Add(FiberHitAna* hit);
    void   SetZ(double _z){i_z = _z;};
    void   SetPos(double _pos){i_pos = _pos;};
    void   SetPosL(double _pos){i_posl = _pos;};
    void   SetdE(double _dE){i_dE = _dE;};
    void   SetUsed(){   i_used = true;};
    void   SetNoUsed(){ i_used = false;};
    void   Print();

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
    double i_pos     = -9999.;
    double i_pos_org = -9999.;
    double i_ang = -9999.;
    double t_tot = -9999.;
    double t_time = -9999.;
    double i_cl_fiber = -9999.;
    double i_z     = -9999.;
    double i_z_org = -9999.;
    double i_res   = -9999.;
    double i_posl  = -9999.;
    double i_dE    = -9999.;
    int  i_did = -1;
    int  i_ud = 0; // 0:upstream 1:downstream
    bool i_valid = false;
    bool i_used  = false;
    int i_simtrackid = -1;
    double GetAngFiber(ParaManager *par);
    double GetPosFiber(ParaManager *par);
    double GetZFiber(  ParaManager *par);
    bool GetValidFiber(ParaManager *par);
};


#endif
