#include "THypernucleus.hh"

using namespace std;

ClassImp(THypernucleus)
THypernucleus::THypernucleus()
{
  Pattern=-1;
  
  //Mother:
  PDG=-1;
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
  Id_Fragment=-1;
  MomE_Fragment.SetXYZT(0.,0.,0.,0.);
  Chi2ndf_Fragment=-1.;
  NDF_Fragment=-999;
  Pvalue_Fragment=-1.;
  Angle_MotherFragment=-1.;
  Fragment_IsFromHyp=-1;

  Id_Pion=-1;
  MomE_Pion.SetXYZT(0.,0.,0.,0.);
  Chi2ndf_Pion=-1.;
  NDF_Pion=-999;
  Pvalue_Pion=-1.;
  Angle_MotherPion=-1.;
  NHitsMDC_Pion=-1;
  NHitsMinifiber_Pion=-1;
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
  PDG=H.PDG;
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
  Id_Fragment=H.Id_Fragment;
  MomE_Fragment=H.MomE_Fragment;
  Chi2ndf_Fragment=H.Chi2ndf_Fragment;
  NDF_Fragment=H.NDF_Fragment;
  Pvalue_Fragment=H.Pvalue_Fragment;
  Angle_MotherFragment=H.Angle_MotherFragment;
  Fragment_IsFromHyp=H.Fragment_IsFromHyp;

  Id_Pion=H.Id_Pion;
  MomE_Pion=H.MomE_Pion;
  Chi2ndf_Pion=H.Chi2ndf_Pion;
  NDF_Pion=H.NDF_Pion;
  Pvalue_Pion=H.Pvalue_Pion;
  Angle_MotherPion=H.Angle_MotherPion;
  NHitsMDC_Pion=H.NHitsMDC_Pion;
  NHitsMinifiber_Pion=H.NHitsMinifiber_Pion;
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
  PDG=-1;
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
  Id_Fragment=-1;
  MomE_Fragment.SetXYZT(0.,0.,0.,0.);
  Chi2ndf_Fragment=-1.;
  NDF_Fragment=-999;
  Pvalue_Fragment=-1.;
  Angle_MotherFragment=-1.;
  Fragment_IsFromHyp=-1;

  Id_Pion=-1;
  MomE_Pion.SetXYZT(0.,0.,0.,0.);
  Chi2ndf_Pion=-1.;
  NDF_Pion=-999;
  Pvalue_Pion=-1.;
  Angle_MotherPion=-1.;
  NHitsMDC_Pion=-1;
  NHitsMinifiber_Pion=-1;
  N_Pion=-1;
  Pion_IsFromHyp=-1;

  Dist_Daughters=-1.;
  ArmPod_Qt=-1.;
  ArmPod_Alfa=-1.;
}

