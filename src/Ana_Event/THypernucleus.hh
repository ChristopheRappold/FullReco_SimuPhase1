#ifndef HYP_h
#define HYP_h

#include "TVector3.h"
#include <string>
#include <vector>
#include "TObject.h"
#include "TLorentzVector.h"

//using namespace std;

class THypernucleus : public TObject{

public :
  
  Int_t Pattern; /// 1 = Simulation / 2 = KFParticle_real / 3 = KFParticle_cut / 4 = KFParticle / 5 = KFParticle_Mass / 6 = LorentzVector
  
  //Mother:
  Int_t PDG;
  Int_t N_Mother;
  Double32_t Chi2ndf;
  TLorentzVector MomE;
  TVector3 PrimVtx;
  TVector3 DecayVtx;
  TVector3 Dist_RealReconsVtx;
  Double32_t Dist_MotherPrimVtx;
  Double32_t Angle_MotherPrimVtx;
  Double32_t InvMass;
  Double32_t ErrInvMass;
  Int_t ErrGetMass;
  Double32_t LifeTime;
  Double32_t ErrLifeTime;
  Int_t ErrGetLifeTime;

  //Daughters:
  Int_t Id_Fragment;
  TLorentzVector MomE_Fragment;
  Double32_t Chi2ndf_Fragment;
  Int_t NDF_Fragment;
  Double32_t Pvalue_Fragment;
  Double32_t Angle_MotherFragment;
  Int_t Fragment_IsFromHyp; // 0-> No; 1-> Yes
  
  Int_t Id_Pion;
  TLorentzVector MomE_Pion;
  Double32_t Chi2ndf_Pion;
  Int_t NDF_Pion;
  Double32_t Pvalue_Pion;
  Double32_t Angle_MotherPion;
  Int_t NHitsMDC_Pion;
  Int_t NHitsMinifiber_Pion;
  Int_t N_Pion;
  Int_t Pion_IsFromHyp; // 0-> No; 1-> Yes

  Double32_t Dist_Daughters;
  Double32_t ArmPod_Qt;
  Double32_t ArmPod_Alfa;

  THypernucleus();
  THypernucleus(const THypernucleus& H);
  ~THypernucleus();

  virtual void Clear(Option_t *option ="");

  ClassDef(THypernucleus,2)
    
};

#endif
