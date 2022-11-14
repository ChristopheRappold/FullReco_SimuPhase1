#ifndef PSFE_HIT_ANA_HH
#define PSFE_HIT_ANA_HH
#include "../EventWASAUnpack/WASAUnpackBranch.hh"
#include "ParaManager.hh"

class PSFEHitAna
{
  public:
    PSFEHitAna(int seg, int t, int q, ParaManager *par);
    ~PSFEHitAna();
    int GetSeg(){return _seg;};
    int GetT(){return _t;};
    double GetRmin(){  return _rmin;};
    double GetRmax(){  return _rmax;};
    double GetPhi(){return _phi;};
    double GetZ(){  return _z;};
    double GetTime(){  return _time;};
    void Print();

  private:
    int _seg   = -1;
    int _t     = -99999.;
    int _q     = -99999.;
    double _phi  = -99999.;
    double _rmin = -99999.;
    double _rmax = -99999.;
    double _z    = -99999.;
    double _time = -99999.;
    bool _flag_fiber = false;
    void SetZRPhi(ParaManager *par);
};


#endif
