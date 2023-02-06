#ifndef OPTICS_MOM_HH
#define OPTICS_MOM_HH
#include "../EventWASAUnpack/WASAUnpackBranch.hh"
#include "ParaManager.hh"
#include "FiberTrackAna.hh"
#include "MWDCTracking.hh"
#include "S4SciHitAna.hh"
#include "TVector3.h"
#include "FiberTrackAna.hh"

class OpticsMom
{
  public:
    OpticsMom(FiberTrackAna *fiber, MWDCTracking *mwdc, S4SciHitAna *s4scin, ParaManager *par);
    ~OpticsMom();
    std::map<std::string, std::vector<float>> par_optics;
    float  GetX4(){     return X4;};
    float  GetA4(){     return A4;};
    float  GetY4(){     return Y4;};
    float  GetB4(){     return B4;};
    float  GetX2(){     return X2;};
    float  GetA2(){     return A2;};
    float  GetY2(){     return Y2;};
    float  GetB2(){     return B2;};
    float  GetZ2(){     return Z2;};
    float  GetY2Rec();
    int    GetPID(){    return PID;};
    double GetTOFCorr(){return tof_corr;};
    double GetMass(){   return mass;};
    double GetMom() {   return momentum;}; //pz calculate from Brho
    TVector3 GetMomV() {return v_mom;};
    void   GetPref();
    FiberTrackAna* GetTrack(){return track;};
    void   CalDelta();
    void   CalMom();
    void   CalVec();
    float  GetA2Rec();
    float  GetB2Rec();
    void   CorrectTOF();
    void   SetBest(){ flag_best = true;};
    bool   IsBest(){  return flag_best;};
    void   Print();
  private:
    float X4 = -9999.;
    float A4 = -9999.;
    float Y4 = -9999.;
    float B4 = -9999.;
    float chi2_s4 = -9999.;
    float X2 = -9999.;
    float A2 = -9999.;
    float Y2 = -9999.;
    float B2 = -9999.;
    float Z2 = -9999.;
    float chi2_dft = -9999.;
    int PID = INT_MIN;
    double tof_sc3141 = -9999.;
    double tof_corr= -9999. ;
    float mass     = -9999. ;
    float momentum = -9999. ;
    float delta    = -9999. ;
    double Pref = -9999.;
    TVector3 v_mom;
    bool flag_best = false;
    FiberTrackAna *track = nullptr;
};

#endif
