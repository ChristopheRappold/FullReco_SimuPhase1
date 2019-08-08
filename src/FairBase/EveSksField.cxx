#include "EveSksField.h"

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


EveSksField::EveSksField(const TString& namefield,bool Telsa,double fac,bool secondMagnet,double SBz,bool Solenoid):TEveMagField()
{

  if(Telsa)
    telsa = 0.1;
  else
    telsa = 0.1;

  if(Solenoid==false)
    {
      FieldMap = new FrsSksHypFieldMapFull(namefield);//AladinFieldMap.root");
      //dynamic_cast<FrsSksHypFieldMapFull*>(FieldMap)->SetPosition((1894.014-100)*0.1,0,(414.+40+20+400.+419.89)*0.1,90.*TMath::Pi()/180.,(90.+25+10)*TMath::Pi()/180.);//->SetPosition(0.,0.,0.,0.); // Geometry already in the magnet center
      dynamic_cast<FrsSksHypFieldMapFull*>(FieldMap)->SetPositionFromGeoManager("TMsystem");
      dynamic_cast<FrsSksHypFieldMapFull*>(FieldMap)->Init();
      if(secondMagnet)
	{
	  dynamic_cast<FrsSksHypFieldMapFull*>(FieldMap)->SetPositionSecondFromGeoManager("TMsystem3","SecondMag");
	  dynamic_cast<FrsSksHypFieldMapFull*>(FieldMap)->InitSecond(0.,0.,SBz); // kG
	}
      dynamic_cast<FrsSksHypFieldMapFull*>(FieldMap)->Print();
      
    }
  else
    {
      FieldMap = new FrsSolenoidHypField();
      dynamic_cast<FrsSolenoidHypField*>(FieldMap)->SetPositionFromGeoManager("TransSetDet","FieldSolenoid");
      dynamic_cast<FrsSolenoidHypField*>(FieldMap)->SetField(0.,fac,SBz*0.1);
      dynamic_cast<FrsSolenoidHypField*>(FieldMap)->Print();
    }

  factor = fac;
  if(TMath::Abs(fac-1.)>1e-6)
    normalized = true;
  else
    normalized = false;

  if(Solenoid)
    {
      normalized = false;
      factor = 1.;
    }
}

EveSksField::~EveSksField()
{
  if(FieldMap)
    {
      delete FieldMap;
      FieldMap= 0;
    }
}

TEveVector EveSksField::GetField(Float_t x, Float_t y, Float_t z) const
{
  Double_t B[3];
  Double_t pos[3];
  
  pos[0]=x;
  pos[1]=y;
  pos[2]=z;
  
  FieldMap->Field(pos,B);
  
  TEveVector temp_B(B[0]*telsa,B[1]*telsa,B[2]*telsa);
  
  // if(scale)
  //   {
  //     double scaling_factor=F_scale(z,p0,p1,p0b,p1b);
  //     temp_B*=scaling_factor;
  //   }

   if(normalized)
     temp_B*=factor;

  // //temp_B*=0.7447/0.75;
  
  return temp_B;
}
