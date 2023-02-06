#ifndef MWDC_HIT_ANA_HH
#define MWDC_HIT_ANA_HH
#include "../EventWASAUnpack/WASAUnpackBranch.hh"
#include "ParaManager.hh"
#include <climits>

class MWDCHitAna
{
  public:
    MWDCHitAna(MWDCHit a, ParaManager *par, int tref);
    ~MWDCHitAna();
    /// define functions
    int GetCtdc(){return i_ctdc;};
    int GetCh(){return i_ch;};
    int GetPlane(){return i_plane;};
    int GetWire(){return i_wire;};
    int GetLT(){return t_leading;};
    int GetTT(){return t_trailing;};
    int GetTOT(){return tot;};
    double GetDt(){return i_dt;};
    double GetDl(){return i_dl;};

  private:
    int i_ctdc = -1;
    int i_ch = -1;
    int i_plane = -1;
    int i_wire = -1;
    int t_leading = INT_MIN;
    int t_trailing = INT_MIN;
    int tot = INT_MIN;
    double i_dt = -9999.;
    double i_dl = -9999.;
    bool i_valid = false;
    int i_did = -1;
    void SetDriftTime(ParaManager *par);
    void SetDriftLength(ParaManager *par);
};
#endif
