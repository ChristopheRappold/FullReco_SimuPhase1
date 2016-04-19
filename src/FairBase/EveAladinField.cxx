#include "EveAladinField.h"

double F_scale(double x,double p0,double p1,double p0b,double p1b)
{
  if(x<60. || x>200.)
    return 1.;

  double v_fit = TMath::Exp( p0+p1*x);
  double v_new = TMath::Exp( p0b+p1b*x);
  if(TMath::Abs(v_fit)<1e-9)
    return 1.;
  else
    return v_new/v_fit;
}


EveAladinField::EveAladinField(bool normalized_,double fac,bool Telsa,bool DoFieldScaling,double p1_scaling):TEveMagField(),normalized(normalized_),scale(DoFieldScaling),factor(fac)
{

  if(Telsa)
    telsa = 0.01;
  else
    telsa = 1.;

  p0 = 1.67525;
  p1 = -3.82014e-2;

  p1b = p1_scaling;
  p0b=p0+(p1-p1b)*60.;

  double Aladin_Angle = -5.6;    //degree
  FieldMap = new HypFieldMapFull("../input/AladinFieldMapFull.root");//AladinFieldMap.root");
  dynamic_cast<HypFieldMapFull*>(FieldMap)->SetPosition(0.,0.,0.,Aladin_Angle*TMath::Pi()/180.);//->SetPosition(0.,0.,0.,0.); // Geometry already in the magnet center
  dynamic_cast<HypFieldMapFull*>(FieldMap)->Init();
  dynamic_cast<HypFieldMapFull*>(FieldMap)->Print();


}

EveAladinField::EveAladinField(TString namefield,bool normalized_,double fac,bool Telsa,bool DoFieldScaling,double p1_scaling):TEveMagField(),normalized(normalized_),scale(DoFieldScaling),factor(fac)
{

  if(Telsa)
    telsa = 0.01;
  else
    telsa = .1;

  p0 = 1.67525;
  p1 = -3.82014e-2;

  p1b = p1_scaling;
  p0b=p0+(p1-p1b)*60.;

  double Aladin_Angle = -5.6;    //degree
  FieldMap = new HypFieldMapFull(namefield.Data());//AladinFieldMap.root");
  dynamic_cast<HypFieldMapFull*>(FieldMap)->SetPosition(0.,0.,0.,Aladin_Angle*TMath::Pi()/180.);//->SetPosition(0.,0.,0.,0.); // Geometry already in the magnet center
  dynamic_cast<HypFieldMapFull*>(FieldMap)->Init();
  dynamic_cast<HypFieldMapFull*>(FieldMap)->Print();


}



EveAladinField::~EveAladinField()
{
  if(FieldMap)
    {
      delete FieldMap;
      FieldMap= 0;
    }
}

TEveVector EveAladinField::GetField(Float_t x, Float_t y, Float_t z) const
{
  Double_t B[3];
  Double_t pos[3];
  
  pos[0]=x;
  pos[1]=y;
  pos[2]=z;
  
  FieldMap->Field(pos,B);
  
  TEveVector temp_B(B[0]*telsa,B[1]*telsa,B[2]*telsa);
  
  if(scale)
    {
      double scaling_factor=F_scale(z,p0,p1,p0b,p1b);
      temp_B*=scaling_factor;
    }

  if(normalized)
    temp_B*=factor/0.75;

  //temp_B*=0.7447/0.75;
  
  return temp_B;
}
