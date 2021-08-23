#ifndef EVEWASAFIELD_H
#define EVEWASAFIELD_H 1


#include "WasaSolenoidFieldMap.h"
#include "TEveTrackPropagator.h"

#include <memory>

class EveWasaField : public TEveMagField
{
  std::unique_ptr<WasaSolenoidFieldMap> FieldMap;

public:
  EveWasaField(const TString& nameField, double maxF, double signD);
  ~EveWasaField();

  using   TEveMagField::GetField;
  
  TEveVector GetField(Float_t x, Float_t y, Float_t z) const;

};


#endif
