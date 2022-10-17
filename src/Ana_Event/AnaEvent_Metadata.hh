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

  uint Unpack_RunNumber = 9999;
  uint Unpack_FirstFile = 9999;
  uint Unpack_LastFile = 9999;
  uint Unpack_EventTotal = 0;
  uint Unpack_EventUnpacked = 0;

  std::string Unpack_SetupFiber = "";
  std::string Unpack_ChannelMapT0 = "";
  std::string Unpack_ChannelMapPSFE = "";
  std::string Unpack_ChannelMapPSBE = "";
  std::string Unpack_SetupPSB = "";
  std::string Unpack_DtDxTableMWDC = "";
  std::string Unpack_CellOffsetS4WFD1 = "";
  std::string Unpack_CellOffsetS2WFD1 = "";
  std::string Unpack_CellOffsetS2WFD2 = "";
  std::string Unpack_CellOffsetS2WFD3 = "";
  std::string Unpack_CellOffsetS2WFD4 = "";
  std::string Unpack_CellOffsetS2WFD5 = "";
  std::string Unpack_ChannelMapMDC = "";
  std::string Unpack_PhysicalMapMDC = "";
  std::string Unpack_DriftParamMDC = "";

  void ShowMeta() const;

  ClassDef(AnaEvent_Metadata,2)

};



#endif
