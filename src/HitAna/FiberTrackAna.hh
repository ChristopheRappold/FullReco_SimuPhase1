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
    FiberTrackAna(std::vector<FiberHitAna*> &cont, ParaManager *par);
    FiberTrackAna(FiberHitXUV *hit_xuv, PSBHitAna *hit_psb, ParaManager *par);
    FiberTrackAna(double x_tgt, double y_tgt, PSBHitAna *hit_psb, ParaManager *par);
    ~FiberTrackAna();
    double GetX(){     return _x   ;};
    double GetY(){     return _y   ;};
    double GetXtgt(){  return _xtgt;};
    double GetYtgt(){  return _ytgt;};
    double GetXdet(){  return _xdet;};
    double GetYdet(){  return _ydet;};
    double GetXmft(){  return _xmft;};
    double GetYmft(){  return _ymft;};
    double GetA(){     return _a   ;};
    double GetB(){     return _b   ;};
    double GetChi2(){  return _chi2;};
    int    GetNlayer(){return _nlayer;};
    double GetTOT(){   return _tot;};
    double GetTime(){  return _time;};
    std::vector<FiberHitXUV*> GetContXUV(){ return _cont_xuv;};
    std::vector<FiberHitAna*> GetContHit(){ return _cont_hit;};
    std::map<int, double>     GetResidual(){return _cont_res;};
    std::map<int, std::map<double, double>> GetPosRes(){return _cont_pos_res;};
    std::map<int, std::vector<double>> GetResidualChPos(){return _cont_res_ch_pos;};
    std::map<int, int> GetDidFib(){return _cont_did_fib;};
    void   CorrectMFT(ParaManager *par);
    void   CorrectMFTCombi(ParaManager *par);
    //void   CorrectMFTXY(ParaManager *par);
    void   SetPosL();
    void   SetNoValid(){_valid = false;}
    bool   IsValid(){return _valid;}
    void   SortContHit();
    FiberHitAna* GetHit(int i){ return _cont_hit[i];};
    bool   IsInclusive(FiberTrackAna* track);
    void   SetBest(){ _flag_best = true;};
    bool   IsBest(){ return _flag_best;};
    void   Print();

    bool   IsFlagXUV(){return   _flag_xuv; };
    bool   IsFlagCombi(){return _flag_combi; };
    bool   IsFlagPair(){return  _flag_pair; };
    void   SetFlagCombi(){ _flag_combi = true;};
    void   SetFlagPair(){  _flag_pair  = true;};

    void   SetFlagPSB(){ _flag_psb = true; };
    void   DelFlagPSB(){ _flag_psb = false; };
    bool   IsFlagPSB(){return _flag_psb; };
    void   SetPSBHit(PSBHitAna *hit){ _psb_hit = hit; };
    PSBHitAna* GetPSBHit(){return _psb_hit; };
    void   SetSegPSB(int seg){ _seg_psb = seg; };
    int    GetSegPSB(){return _seg_psb; };
    void   SetPSBDifZ(  double dif_z){ _psb_dif_z = dif_z; };
    void   SetPSBDifPhi(double dif_p){ _psb_dif_p = dif_p; };
    double GetPSBDifZ(){   return _psb_dif_z;};
    double GetPSBDifPhi(){ return _psb_dif_p;};


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
    double _xdet = -9999.;
    double _ydet = -9999.;
    double _xmft = -9999.;
    double _ymft = -9999.;
    double _a    = -9999.;
    double _b    = -9999.;
    double _z_mean = -9999.;
    double _chi2 = -9999.;
    double _tot  = -9999.;
    double _time = -9999.;
    int    _nlayer     = -1;
    bool   _flag_xuv   = false;
    bool   _flag_combi = false;
    bool   _flag_pair  = false;
    bool   _flag_psb   = false;
    int    _seg_psb    = -1;
    double _psb_dif_z  = -9999.;
    double _psb_dif_p  = -9999.;
    bool   _flag_psfe  = false;
    int    _seg_psfe   = -1;
    bool   _flag_psbe  = false;
    int    _seg_psbe   = -1;
    bool   _valid = true;
    bool   _flag_best  = false;
    PSBHitAna*  _psb_hit  = nullptr;
    PSFEHitAna* _psfe_hit = nullptr;
    std::vector<FiberHitXUV*> _cont_xuv;
    std::vector<FiberHitAna*> _cont_hit;
    std::map<int, double>     _cont_res;
    std::map<int, std::vector<double>>  _cont_res_ch_pos;
    std::map<int, std::map<double, double>> _cont_pos_res;
    std::map<int, int> _cont_did_fib;
    void Tracking(std::vector<FiberHitXUV*> &cont, ParaManager *par);
    void Tracking(std::vector<FiberHitAna*> &cont, ParaManager *par);
    void Tracking(FiberHitXUV *hit_xuv, PSBHitAna *hit_psb, ParaManager *par);
    void Tracking(double x_tgt, double y_tgt, PSBHitAna *hit_psb, ParaManager *par);
    bool TrackFitting(int nlayer, double* w, double* z, double* angle, double* s);
    bool GaussJordan( double **a, int n, double *b, int *indxr, int *indxc, int *ipiv );
    void SetResidual();
    void SetTOT();
    void SetTime();
};


#endif
