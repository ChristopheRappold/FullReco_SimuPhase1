#ifndef  MWDCTRACKING_HH
#define  MWDCTRACKING_HH
#include <iostream>
#include "TMath.h"
#include "MWDCHitAna.hh"
//#define  MWDC_RESOLUTION_ASSUMED 0.5 // mm

class MWDCTracking
{
  public:
    MWDCTracking(ParaManager *par);
    ~MWDCTracking();
    //----------------
    void StoreHit(int, Float_t, Float_t, Float_t, Float_t,  Float_t, Float_t);
    void StoreHit(MWDCHitAna *a, ParaManager *par); ///add by Enqiang
    int Tracking(void);
    int Tracking_FitCylindrical(void);
    int Tracking_LRfixed(void);
    int Tracking_PairLR(void);
    int ResidualExCalc(void);
    int ResidualExCalc2(void);
    int ResidualExCalc3(void);
    int FindTime_wMaxToT(void);
    //----------------
    Float_t GetX(void);
    Float_t GetY(void);
    Float_t GetA(void);
    Float_t GetB(void);
    void SetX( Float_t _x_track) { x_track = _x_track; };
    void SetY( Float_t _y_track) { y_track = _y_track; };
    void SetA( Float_t _a_track) { a_track = _a_track; };
    void SetB( Float_t _b_track) { b_track = _b_track; };
    Float_t GetXError(void);
    Float_t GetYError(void);
    Float_t GetAError(void);
    Float_t GetBError(void);
    Float_t GetChi2(void);
    Float_t GetResidual(int);
    Float_t GetTime_wMaxToT(int);
    Float_t GetExclusiveResidual(int);
    void SetAssumedResolution(float);
    void SetMaxHitCombination(int);
    void SetMinPlaneEnebaled(int);
    int  GetTrackingStatus(void);
    int  GetNumEnabledPlane(void);
    int  GetNumLRCombination(void);
    int  GetNumHitCombination(void);
    void GetPlaneEnabled(int*);
    void GetHitCombination(int*);
    void GetLRCombination(int*);
    int  GetLR(int);
    int  GetiHit(int);
    Float_t  GetUsedWirepos(int);
    Float_t  GetUsedLength(int);
    //----------------


  private:
    int TrackingInverseMatrix(void);
    int TrackingFitCylindrical(void);
    int NumEnabledPlaneTemp(void);
    int NumLRCombinationTemp(void);
    int NumHitCombinationTemp(void);
    int HitCombinationTempIncrement(void);
    int LRCombinationTempIncrement(void);
    void SetHitCombinationTemp0(void);
    void SetLRCombinationTemp0(void);
    void CopyPlaneEnabled(int*, int*);
    void CopyHitCombination(int*, int*);
    void CopyLRCombination(int*, int*);
    int SetPairTemp(void);
    int NumPairCombinationTemp(void);
    int SetPairCombinationTemp0(void);
    int SetHitLRfromPair(void);
    int PairCombinationTempIncrement(void);
    int CheckSameType(int,int);
    //---------------
    //
    Float_t mwdc_tracking_assumption_resoltion;
    int numhit[16];
    Float_t wirepos[16][64];
    Float_t length[16][64];
    Float_t tot[16][64];
    Float_t tot_max[16];
    Float_t time[16][64];
    Float_t theta[16];
    Float_t zplane[16];
    int max_hit_combination;
    int min_plane_enabled;
    int min_uv_plane_enabled;
    //-------------
    Float_t chi2_temp;
    double x_temp, a_temp, y_temp, b_temp;
    double xerr_temp, aerr_temp, yerr_temp, berr_temp;
    Float_t residual_temp[16];
    int pair_enabled_temp[8];
    int plane_enabled_temp[16], i_hit_used_temp[16], LR_used_temp[16];
    bool LR_fixed[16];
    //-------------
    Float_t x_track, a_track, y_track, b_track;
    Float_t xerr_track, aerr_track, yerr_track, berr_track;
    Int_t time_wMaxToT[16];
    Float_t chi2_track;
    Float_t residual_track[16];
    Float_t residual_track_ex[16];
    int tracking_status;
    int plane_enabled_track[16], i_hit_used_track[16], LR_used_track[16];
    //--------------
    int i_hit1_pair[8][32];
    int i_hit2_pair[8][32];
    int i_lr1_pair[8][32];
    int i_lr2_pair[8][32];
    int num_pair[8];
    int i_pair_used_temp[8];
    int max_pair_combination;
    //
};

#endif // MWDCTRACKING_HH
