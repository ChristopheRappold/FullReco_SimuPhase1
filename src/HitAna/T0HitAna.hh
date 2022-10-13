#ifndef T0_HIT_ANA_HH
#define T0_HIT_ANA_HH
#include "../EventWASAUnpack/WASAUnpackBranch.hh"
#include "ParaManager.hh"

class T0HitAna
{
  public:
    T0HitAna(int seg, int t_u, int t_d, int q_u, int q_d, ParaManager *par);
    ~T0HitAna();
    int GetSeg(){return _seg;};
    int GetTU(){return _t_u;};
    int GetTD(){return _t_d;};
    int GetQU(){return _q_u;};
    int GetQD(){return _q_d;};
    double GetZ(){  return _z;};
    double GetTimeU(){  return _time_u;};
    double GetTimeD(){  return _time_d;};
    double GetTime(){   return _time;};
    double GetDT(){     return _time_u - _time_d;};
    void Print();

  private:
    int _seg   = -1;
    int _t_u   = -99999.;
    int _t_d   = -99999.;
    int _q_u   = -99999.;
    int _q_d   = -99999.;
    double _z   = -99999.;
    double _time_u = -99999.;
    double _time_d = -99999.;
    double _time = -99999.;
    void SetZ(ParaManager *par);
};


#endif
