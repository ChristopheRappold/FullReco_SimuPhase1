#ifndef EVESKSFIELD_H
#define EVESKSFIELD_H 1


#include "FrsSksHypFieldMapFull.h"
#include "FrsSolenoidHypField.h"
#include "TEveTrackPropagator.h"


class EveSksField : public TEveMagField
{
  FairField* FieldMap;  

  double telsa;
  double factor;
  bool normalized;

public:
  EveSksField(const TString& n,bool Telsa=false,double fac =1.,bool SecondM = false,double Sbz = 10.,bool Solenoid=false);
  ~EveSksField();

  using   TEveMagField::GetField;
  
  TEveVector GetField(Float_t x, Float_t y, Float_t z) const;

};


#endif
