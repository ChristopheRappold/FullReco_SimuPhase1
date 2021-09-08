#include "EveWasaField.h"
#include "WasaSolenoidFieldMap.h"

// double F_scale(double x,double p0,double p1,double p0b,double p1b)
// {
//   if(x<60. || x>200.)
//     return 1.;

//   double v_fit = TMath::Exp( p0+p1*x);
//   double v_new = TMath::Exp( p0b+p1b*x);
//   if(TMath::Abs(v_fit)<1e-9)
//     return 1.;
//   else
//     return v_new/v_fit;
// }


EveWasaField::EveWasaField(const TString& namefield, double maxField, double signD):TEveMagField()
{
  FieldMap = std::make_unique<WasaSolenoidFieldMap>("WasaSolEveFieldMap","WasaSolEveFieldMap",namefield,maxField,signD);
  FieldMap->Init();
}

EveWasaField::~EveWasaField()
{ }

TEveVector EveWasaField::GetField(Float_t x, Float_t y, Float_t z) const
{
  Double_t B[3];
  Double_t pos[3];
  
  pos[0]=x;
  pos[1]=y;
  pos[2]=z;
  
  FieldMap->GetBxyz(pos,B);
  
  TEveVector temp_B(B[0]*0.1,B[1]*0.1,B[2]*0.1);  // in tesla

  return temp_B;
}
