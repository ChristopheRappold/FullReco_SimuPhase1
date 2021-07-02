#include "THypernucleus.hh"

using namespace std;

ClassImp(THypernucleus)
THypernucleus::THypernucleus()
{
  Pattern=-1;
  
  //Mother:
  N_Mother=-1;
  Chi2ndf=-1.;
  MomE.SetPxPyPzE(0.,0.,0.,0.);
  PrimVtx.SetXYZ(0.,0.,0.);
  DecayVtx.SetXYZ(0.,0.,0.);
  Dist_RealReconsVtx.SetXYZ(-100.,-100.,-100.);
  Dist_MotherPrimVtx=-1.;
  Angle_MotherPrimVtx=-1.;
  InvMass=-1.;
  ErrInvMass=-1.;
  ErrGetMass=-1;
  LifeTime=-1.;
  ErrLifeTime=-1.;
  ErrGetLifeTime=-1;


  //Daughters:
  MomE_Fragment.SetXYZT(0.,0.,0.,0.);
  Angle_MotherFragment=-1.;
  MomE_Pion.SetXYZT(0.,0.,0.,0.);
  Chi2ndf_Pion=-1.;
  Angle_MotherPion=-1.;
  N_Pion=-1;
  Pion_IsFromHyp=-1;
  Dist_Daughters=-1.;
  ArmPod_Qt=-1.;
  ArmPod_Alfa=-1.;
  }

THypernucleus::THypernucleus(const THypernucleus& H)
{
  Pattern=H.Pattern;
  
  //Mother:
  N_Mother=H.N_Mother;
  Chi2ndf=H.Chi2ndf;
  MomE=H.MomE;
  PrimVtx=H.PrimVtx;
  DecayVtx=H.DecayVtx;
  Dist_RealReconsVtx=H.Dist_RealReconsVtx;
  Dist_MotherPrimVtx=H.Dist_MotherPrimVtx;
  Angle_MotherPrimVtx=H.Angle_MotherPrimVtx;
  InvMass=H.InvMass;
  ErrInvMass=H.ErrInvMass;
  ErrGetMass=H.ErrGetMass;
  LifeTime=H.LifeTime;
  ErrLifeTime=H.ErrLifeTime;
  ErrGetLifeTime=H.ErrGetLifeTime;

  //Daughters:
  MomE_Fragment=H.MomE_Fragment;
  Angle_MotherFragment=H.Angle_MotherFragment;
  MomE_Pion=H.MomE_Pion;
  Chi2ndf_Pion=H.Chi2ndf_Pion;
  Angle_MotherPion=H.Angle_MotherPion;
  N_Pion=H.N_Pion;
  Pion_IsFromHyp=H.Pion_IsFromHyp;
  Dist_Daughters=H.Dist_Daughters;
  ArmPod_Qt=H.ArmPod_Qt;
  ArmPod_Alfa=H.ArmPod_Alfa;
}

//THypernucleus& THypernucleus::operator=(const THypernucleus& H)
//{
//   THypernucleus temp(H);
//   return temp;
//}


THypernucleus::~THypernucleus()
{
}

void THypernucleus::Clear(Option_t *option)
{
  Pattern=-1;
  
  //Mother:
  N_Mother=-1;
  Chi2ndf=-1.;
  MomE.SetPxPyPzE(0.,0.,0.,0.);
  PrimVtx.SetXYZ(0.,0.,0.);
  DecayVtx.SetXYZ(0.,0.,0.);
  Dist_RealReconsVtx.SetXYZ(-100.,-100.,-100.);
  Dist_MotherPrimVtx=-1.;
  Angle_MotherPrimVtx=-1.;
  InvMass=-1.;
  ErrInvMass=-1.;
  ErrGetMass=-1;
  LifeTime=-1.;
  ErrLifeTime=-1.;
  ErrGetLifeTime=-1;

  //Daughters:
  MomE_Fragment.SetXYZT(0.,0.,0.,0.);
  Angle_MotherFragment=-1.;
  MomE_Pion.SetXYZT(0.,0.,0.,0.);
  Chi2ndf_Pion=-1.;
  Angle_MotherPion=-1.;
  N_Pion=-1;
  Pion_IsFromHyp=-1;
  Dist_Daughters=-1.;
  ArmPod_Qt=-1.;
  ArmPod_Alfa=-1.;
}

