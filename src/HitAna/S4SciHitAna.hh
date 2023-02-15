#ifndef S4SCI_HIT_ANA_HH
#define S4SCI_HIT_ANA_HH
#include "../EventWASAUnpack/WASAUnpackBranch.hh"
#include "ParaManager.hh"
#include "../FullRecoConfig.hh"

class S4SciHitAna
{
  public:
    S4SciHitAna(S4TQ *a, ParaManager *par, const std::string& StudyCase); /// modify later
    ~S4SciHitAna();

    float GetTOF_sc3141(){return tof_sc3141;};
    float GetTOF_sc4143(){return tof_sc4143;};
    float GetdE_sc31(){return de_sc31;};
    float GetdE_sc41(){return de_sc41;};
    float GetdE_sc42_low(){return de_sc42_low;};
    float GetdE_sc42_high(){return de_sc42_high;};
    int   GetPID(){ return PID_residual;};
    float GetMass(){return mass_residual;};
    void  ResIdentify(const std::string& StudyCase);
    void  Print();

  private:
    int tdc_sc31[2];
    int tdc_sc41[3][2];
    int tdc_sc42[2];
    int tdc_sc43[4][2];
    int tdc_trig[3];
    int tdc_s2tref[3];
    int qdc_sc31[2];
    int qdc_sc41[3][2];
    int qdc_sc42[2];
    int qdc_sc42_lowgain[2];
    int qdc_sc43[4][2];
    float t_sc31 = -999.;
    float t_sc41 = -999.;
    float t_sc42 = -999.;
    float t_sc43 = -999.;
    float de_sc31 = -999.;
    float de_sc41 = -999.;
    float de_sc42_low = -999.;
    float de_sc42_high = -999.;
    float de_sc43 = -999.;
    float tof_sc3141 = -999.;
    float tof_sc4143 = -999.;
    float t_trig = -999.;
    float t_s2tref = -999.;
    int PID_residual = -999;
    float mass_residual = 0.;
};
#endif
