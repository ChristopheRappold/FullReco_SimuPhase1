#include "THyphiTrack.hh"

using namespace std;

ClassImp(THyphiTrack)
  
THyphiTrack::THyphiTrack():type(""),MC_status(0),Chi2(-1.),Chi2_X(-1.),Chi2_Y(-1.),Mass(-1.),pdgcode(0),MomMass(0.,0.,0.,0.),Mom(0.,0.,0.),Charge(0),BarId(0),Beta(0),HitTR1(-999,-999,-999),HitTR2(-999,-999,-999),HitDC2(-999,-999,-999),NumDC2(-1,-1,-1),Pval2(-1),TofsBar(-1)
{
}

THyphiTrack::THyphiTrack(const THyphiTrack& H)
{
  type=H.type;
  MC_status=H.MC_status;
  Chi2=H.Chi2;
  Chi2_X=H.Chi2_X;
  Chi2_Y=H.Chi2_Y;
  Mass=H.Mass;
  pdgcode=H.pdgcode;
  MomMass=H.MomMass;
  Mom=H.Mom;
  
  Charge=H.Charge;
  BarId =H.BarId;
  dE = H.dE;
  Beta = H.Beta;
  RefPoint = H.RefPoint;
  HitTR1 = H.HitTR1;
  HitTR2 = H.HitTR2;
  HitDC2 = H.HitDC2;
  NumDC2 = H.NumDC2;
  HitTOF = H.HitTOF;
  Pval2  = H.Pval2;
  TofsBar= H.TofsBar;
}

//THyphiTrack& THyphiTrack::operator=(const THyphiTrack& H)
//{
//   THyphiTrack temp(H);
//   return temp;
//}
void THyphiTrack::SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec)
{
  type=name;
  MC_status=MC;
  Chi2 = chi2;
  Chi2_X = chi2_x;
  Chi2_Y = chi2_y;
  Mass = mass;
  pdgcode = pdg;
  MomMass.SetXYZM(vec.X(),vec.Y(),vec.Z(),mass);
  Mom = vec;
  

}
void THyphiTrack::SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec,Int_t charge,Int_t barid,Double_t beta,const TVector3& hit_tr1,const TVector3& hit_tr2, const TVector3 hit_dc2,const TVector3 num_dc2,Double_t pval2,Int_t tofsbar )
{

  type=name;
  MC_status=MC;
  Chi2 = chi2;
  Chi2_X = chi2_x;
  Chi2_Y = chi2_y;
  Mass = mass;
  pdgcode = pdg;
  MomMass.SetXYZM(vec.X(),vec.Y(),vec.Z(),mass);
  Mom = vec;

  Charge=charge;
  BarId=barid;
  Beta=beta;
  HitTR1=hit_tr1;
  HitTR2=hit_tr2;
  HitDC2=hit_dc2;
  NumDC2 = num_dc2;
  Pval2 =pval2;
  TofsBar=tofsbar;

}

THyphiTrack::~THyphiTrack()
{
}

void THyphiTrack::Clear(Option_t *option)
{

  MomMass.SetXYZT(0.,0.,0.,0.);
  Mom.SetXYZ(0.,0.,0.);

  RefPoint.SetXYZ(-999.,-999.,-999.);

  HitTR1.SetXYZ(-999.,-999.,-999.);
  HitTR2.SetXYZ(-999.,-999.,-999.);
  HitDC2.SetXYZ(-999.,-999.,-999.);

  NumDC2.SetXYZ(-1,-1,-1);

  HitTOF.SetXYZ(-999.,-999.,-999.);

  type="";
  MC_status=-999;
  Chi2 = -999.;
  Chi2_X = -999.;
  Chi2_Y = -999.;
  Mass = -999.;
  pdgcode = -999;

  Charge=-999;
  BarId=-999;
  dE = -999;
  Beta=-999.;
  Pval2 = -999;
  TofsBar=-1;

}

