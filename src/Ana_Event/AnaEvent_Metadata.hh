#ifndef AnaEvent_Metadata_h
#define AnaEvent_Metadata_h

#include <string>

#include "TObject.h"

class AnaEvent_Metadata : public TObject
{
 public :
  std::string NameIn = "";
  std::string NameOut = "";
  std::string DateOfRun = "";
  std::string Hash = "";
  std::string FirstStep = "";
  std::string FinalStep = "";

  bool G4_simu = true;

  uint NEvent = 0;
  uint StartEvent = 0;
  uint StopEvent = 0;
  uint Nb_Fraction = 0;

  uint Wasa_Side = 1;
  bool Wasa_FieldMap = false;
  double Field_Strength = 1.;
  std::string Wasa_FieldMapName = "";

  void ShowMeta() const;

  ClassDef(AnaEvent_Metadata,2)

};



#endif
