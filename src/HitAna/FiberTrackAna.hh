#ifndef FIBER_TRACK_ANA_HH
#define FIBER_TRACK_ANA_HH
#include "FiberHitXUV.hh"
#include "PSBHitAna.hh"
#include "PSFEHitAna.hh"
#include <map>

class FiberTrackAna
{
  public:
    FiberTrackAna(std::vector<FiberHitXUV*> &cont, ParaManager *par);
    ~FiberTrackAna();
    double GetX(){    return _x   ;};
    double GetY(){    return _y   ;};
    double GetXtgt(){ return _xtgt;};
    double GetYtgt(){ return _ytgt;};
    double GetA(){    return _a   ;};
    double GetB(){    return _b   ;};
    double GetChi2(){ return _chi2;};
    std::vector<FiberHitXUV*> GetContXUV(){ return _cont_xuv;};
    std::map<int, double>     GetResidual(){return _cont_res;};
    std::map<int, std::vector<double>> GetResidualChPos(){return _cont_res_ch_pos;};
    void   Print();
    void   SetFlagPSB(){ _flag_psb = true; };
    bool   IsFlagPSB(){return _flag_psb; };
    void   SetPSBHit(PSBHitAna *hit){ _psb_hit = hit; };
    PSBHitAna* GetPSBHit(){return _psb_hit; };
    void   SetSegPSB(int seg){ _seg_psb = seg; };
    int    GetSegPSB(){return _seg_psb; };
    void   SetFlagPSFE(){ _flag_psfe = true; };
    bool   IsFlagPSFE(){return _flag_psfe; };
    void   SetSegPSFE(int seg){ _seg_psfe = seg; };
    int    GetSegPSFE(){return _seg_psfe; };
    void   SetPSFEHit(PSFEHitAna *hit){ _psfe_hit = hit; };
    PSFEHitAna* GetPSFEHit(){return _psfe_hit; };
    void   SetFlagPSBE(){ _flag_psbe = true; };
    bool   IsFlagPSBE(){return _flag_psbe; };
    void   SetSegPSBE(int seg){ _seg_psbe = seg; };
    int    GetSegPSBE(){return _seg_psbe; };

  private:
    double _x    = -9999.;
    double _y    = -9999.;
    double _xtgt = -9999.;
    double _ytgt = -9999.;
    double _a    = -9999.;
    double _b    = -9999.;
    double _chi2 = -9999.;
    bool   _flag_psb = false;
    int    _seg_psb = -1;
    bool   _flag_psfe = false;
    int    _seg_psfe = -1;
    bool   _flag_psbe = false;
    int    _seg_psbe = -1;
    PSBHitAna*  _psb_hit  = nullptr;
    PSFEHitAna* _psfe_hit = nullptr;
    std::vector<FiberHitXUV*> _cont_xuv;
    std::map<int, double>     _cont_res;
    std::map<int, std::vector<double>>  _cont_res_ch_pos;
    void Tracking(std::vector<FiberHitXUV*> &cont);
    bool TrackFitting(int nlayer, double* w, double* z, double* angle, double* s);
    bool GaussJordan( double **a, int n, double *b, int *indxr, int *indxc, int *ipiv );
    void SetResidual();
};


#endif
