#ifndef EVEALADINFIELD_H
#define EVEALADINFIELD_H 1


#include "HypFieldMap.h"
#include "HypFieldMapFull.h"
#include "TEveTrackPropagator.h"


class EveAladinField : public TEveMagField
{
  FairField* FieldMap;  
  bool normalized;
  bool scale;

  double factor;
  double telsa;
  // exp(p0b + p1b*x) with p0b = p0 + (p1-p1b)*x0 and x0 = 60.
  double p0b;
  double p1b;
  // exp(p0+p1*x) : tail scaling 
  double p0;
  double p1;
  
public:
  EveAladinField(bool normalized_=true,double fac=0.7447,bool Telsa=false,bool DoFieldScaling=false,double p1_scaling=-3.7e-2);
  EveAladinField(const TString& name_field,bool normalized_=true,double fac=0.7447,bool Telsa=false,bool DoFieldScaling=false,double p1_scaling=-3.7e-2);
  ~EveAladinField();

  using   TEveMagField::GetField;
  
  TEveVector GetField(Float_t x, Float_t y, Float_t z) const;

};


#endif
