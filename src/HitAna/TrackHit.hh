#ifndef TRACK_HIT_HH
#define TRACK_HIT_HH
#include "FiberHitAna.hh"
#include "FiberTrackAna.hh"
#include "MDCHitAna.hh"
#include "PSBHitAna.hh"
#include "PSFEHitAna.hh"
#include "ParaManager.hh"

class TrackHit
{
  public:
    TrackHit();
    ~TrackHit();
    void Add(FiberHitAna *hit_fiber);
    void Add(MDCHitAna *hit_mdc);
    void Add(PSBHitAna *hit_psb){   _psb_hit  = hit_psb; };
    void Add(PSFEHitAna *hit_psfe){ _psfe_hit = hit_psfe; };
    void DeleteDupMDC();
    int  GetNum(){     return _fiber_cont.size() + _mdc_cont.size();};
    int  GetNumFiber(){return _fiber_cont.size();};
    int  GetNumMDC(){  return _mdc_cont.size();};
    FiberHitAna*  GetFiberHit(int i);
    MDCHitAna*    GetMDCHit(int i);
    std::vector<FiberHitAna*>  GetFiberHitCont(){return _fiber_cont;};
    std::vector<MDCHitAna*>    GetMDCHitCont(){  return _mdc_cont  ;};
    PSBHitAna*    GetPSBHit(){  return _psb_hit  ;};
    PSFEHitAna*   GetPSFEHit(){ return _psfe_hit  ;};
    void SetTrack(FiberTrackAna *track, ParaManager *par);
    double GetTrackX(){return _x;};
    double GetTrackY(){return _y;};
    double GetTrackZ(){return _z;};
    double GetTrackA(){return _a;};
    double GetTrackB(){return _b;};
    void SetFlagPSB(){ _flag_psb = true; };
    bool IsFlagPSB(){return _flag_psb; };
    void SetFlagPSFE(){ _flag_psfe = true; };
    bool IsFlagPSFE(){return _flag_psfe; };
    void Print();

  private:
    std::vector<FiberHitAna*> _fiber_cont;
    std::vector<MDCHitAna*>   _mdc_cont;
    PSBHitAna*   _psb_hit;
    PSFEHitAna*  _psfe_hit;
    double _x = -9999.;
    double _y = -9999.;
    double _z = -9999.;
    double _a = -9999.;
    double _b = -9999.;
    bool   _flag_psb  = false;
    bool   _flag_psfe = false;
};


#endif
