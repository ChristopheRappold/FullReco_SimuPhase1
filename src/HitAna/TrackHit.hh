#ifndef TRACK_HIT_HH
#define TRACK_HIT_HH
#include "FiberHitAna.hh"
#include "FiberTrackAna.hh"
#include "MDCHitAna.hh"
#include "PSBHitAna.hh"
#include "PSFEHitAna.hh"
#include "ParaManager.hh"
#include "TVector3.h"
#include <unordered_map>

class TrackHit
{
  public:
    TrackHit();
    ~TrackHit();
    void AddFiber(int layer, int fiberID);
    void AddPSB(int psbID);
    void AddPSFE(int psfeID);
    void DeleteDupMDC();
    void SetDidCh();
    int  GetNum(){ return 17 - std::count(_mdchit_bestdif.begin(), _mdchit_bestdif.end(), -1)
                            + 6 - std::count(_fiberhit.begin(), _fiberhit.end(), -1); };
    int  GetNumFiber(){return 6 - std::count(_fiberhit.begin(), _fiberhit.end(), -1);};
    int  GetNumMDC(){  return 17 - std::count(_mdchit_bestdif.begin(), _mdchit_bestdif.end(), -1); };
    std::vector<int> GetFiberHit(){return _fiberhit; };
    std::vector<int> GetMDCHit(){return _mdchit_bestdif; };
    int GetPSBHit(){return _psbhit; };
    int GetPSFEHit(){return _psfehit; };
    void SetTrack(FiberTrackAna *track, double tgt_posz);
    double GetTrackX(){return _x;};
    double GetTrackY(){return _y;};
    double GetTrackZ(){return _z;};
    double GetTrackA(){return _a;};
    double GetTrackB(){return _b;};
    void SetChi2ndf(double chi2ndf){_chi2ndf = chi2ndf; };
    double GetChi2NDF(){return _chi2ndf;};
    void SetFlagPSB(){ _flag_psb = true; };
    bool IsFlagPSB(){return _flag_psb; };
    void SetFlagPSFE(){ _flag_psfe = true; };
    bool IsFlagPSFE(){return _flag_psfe; };
    std::map<int, int> GetDidCh(){return _cont_did_ch;};
    bool IsInclusive(TrackHit* track);
    void SetMDCLayHitCont();
    //void ReplaceMDC(std::map<int, std::pair<TVector3, TVector3> > track_pos_mom, ParaManager *par);
    void SetMDCdif(int lay, int hit, double dif){_mdchit_dif[lay][hit]=dif; };

  private:

    std::vector<std::unordered_map<int,double>> _mdchit_dif;
    std::vector<int> _fiberhit = std::vector<int>(6,-1);
    int _psbhit = -1;
    int _psfehit = -1;
    std::vector<int> _mdchit_bestdif = std::vector<int>(17,-1);
    double _chi2ndf = -9999.;
    double _x = -9999.;
    double _y = -9999.;
    double _z = -9999.;
    double _a = -9999.;
    double _b = -9999.;
    bool   _flag_psb  = false;
    bool   _flag_psfe = false;
    std::map<int, int> _cont_did_ch;
};


#endif
