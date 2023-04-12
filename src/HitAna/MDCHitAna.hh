#ifndef MDC_HIT_ANA_HH
#define MDC_HIT_ANA_HH
#include "../EventWASAUnpack/WASAUnpackBranch.hh"
#include "ParaManager.hh"
#include "limits.h"

class MDCHitAna
{
  public:
    MDCHitAna(MDCHit a, ParaManager *par, int tref, double t_t0);
    ~MDCHitAna();
    int GetCtdc(){   return i_ctdc    ;};
    int GetCh(){     return i_ch      ;};
    int GetLay(){    return i_layer   ;};
    int GetWir(){    return i_wire    ;};
    int GetTL(){     return t_leading ;};
    int GetTT(){     return t_trailing;};
    int GetTot(){    return t_tot;};
    double GetZ(){   return i_z;  };
    double GetR(){   return i_r;  };
    double GetPhi(){ return i_phi;};
    double GetDt(){  return i_dt;};
    double GetDl(){  return i_dl;};
    double GetDpf(){ return i_dpf;};
    double GetDpb(){ return i_dpb;};
    double GetLf(){  return i_lf;};
    double GetLb(){  return i_lb;};
    double GetZf(){  return i_zf;};
    double GetZb(){  return i_zb;};
    double GetWsize(){return i_wsize;};
    double GetRes(){ return i_res;};
    double GetDif(){ return i_dif;};
    void SetDif(double dif){i_dif = dif;};
    int  GetDid(){   return i_did;};
    bool IsValid(){  return i_valid;};
    void Print();

  private:
    int i_ctdc     = -1;
    int i_ch       = -1;
    int i_layer    = -1;
    int i_wire     = -1;
    int t_leading  = INT_MIN;
    int t_trailing = INT_MIN;
    int t_tot      = INT_MIN;
    double i_z   = -9999.;
    double i_r   = -9999.;
    double i_phi = -9999.;
    double i_dt  = -9999.;
    double i_dl  = -9999.;
    double i_dpf = -9999.;
    double i_dpb = -9999.;
    double i_lf  = -9999.;
    double i_lb  = -9999.;
    double i_zf  = -9999.;
    double i_zb  = -9999.;
    double i_wsize = -9999.;
    double i_res  = -9999.;
    double i_dif  = -9999.;
    int  i_did   = -1;
    bool i_valid = false;
    bool GetValidMDC(ParaManager *par);
    void SetDriftTime(ParaManager *par, double t_t0);
    void SetDriftLength(ParaManager *par);
    void SetPhys(ParaManager *par);
};


#endif
