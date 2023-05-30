#include "AnaEvent_Metadata.hh"
#include <iostream>

ClassImp(AnaEvent_Metadata)

void AnaEvent_Metadata::ShowMeta() const
{
  std::cout<<"Metadata :\n";
  std::cout<<"NameIn  :"<<NameIn<<"\n";
  std::cout<<"NameOut :"<<NameOut<<"\n";
  std::cout<<"DateOfRun :"<<DateOfRun<<"\n";
  std::cout<<"Hash :"<<Hash<<"\n";
  std::cout<<"First Step :"<<FirstStep<<"\n";
  std::cout<<"Final Step :"<<FinalStep<<"\n";
  std::cout<<"Geant4 simu :"<<G4_simu<<"\n";
  std::cout<<"NEvent :"<<NEvent<<"\n";
  std::cout<<"StartEvent :"<<StartEvent<<"\n";
  std::cout<<"StopEvent :"<<StopEvent<<"\n";
  std::cout<<"Nb Fraction :"<<Nb_Fraction<<"\n";
  std::cout<<"Wasa Side :"<<Wasa_Side<<"\n";
  std::cout<<"Wasa Field Map :"<<Wasa_FieldMap<<"\n";
  std::cout<<"Wasa Field Map Name :"<<Wasa_FieldMapName<<"\n";
  std::cout<<"Field Strength :"<<Field_Strength<<"\n";


  std::cout<<"Unpack_RunNumber: " << Unpack_RunNumber << "\n";
  std::cout<<"Unpack_FirstFile: " << Unpack_FirstFile << "\n";
  std::cout<<"Unpack_LastFile: " << Unpack_LastFile << "\n";
  std::cout<<"Unpack_EventTotal: " << Unpack_EventTotal << "\n";
  std::cout<<"Unpack_EventUnpacked: " << Unpack_EventUnpacked << "\n";

  std::cout<<"Unpack_SetupFiber: " << Unpack_SetupFiber << "\n";
  std::cout<<"Unpack_ChannelMapT0: " << Unpack_ChannelMapT0 << "\n";
  std::cout<<"Unpack_ChannelMapPSFE: " << Unpack_ChannelMapPSFE << "\n";
  std::cout<<"Unpack_ChannelMapPSBE: " << Unpack_ChannelMapPSBE << "\n";
  std::cout<<"Unpack_SetupPSB: " << Unpack_SetupPSB << "\n";
  std::cout<<"Unpack_DtDxTableMWDC: " << Unpack_DtDxTableMWDC << "\n";
  std::cout<<"Unpack_CellOffsetS4WFD1: " << Unpack_CellOffsetS4WFD1 << "\n";
  std::cout<<"Unpack_CellOffsetS2WFD1: " << Unpack_CellOffsetS2WFD1 << "\n";
  std::cout<<"Unpack_CellOffsetS2WFD2: " << Unpack_CellOffsetS2WFD2 << "\n";
  std::cout<<"Unpack_CellOffsetS2WFD3: " << Unpack_CellOffsetS2WFD3 << "\n";
  std::cout<<"Unpack_CellOffsetS2WFD4: " << Unpack_CellOffsetS2WFD4 << "\n";
  std::cout<<"Unpack_CellOffsetS2WFD5: " << Unpack_CellOffsetS2WFD5 << "\n";
  std::cout<<"Unpack_ChannelMapMDC: " << Unpack_ChannelMapMDC << "\n";
  std::cout<<"Unpack_PhysicalMapMDC: " << Unpack_PhysicalMapMDC << "\n";
  std::cout<<"Unpack_DriftParamMDC: " << Unpack_DriftParamMDC << "\n";
}
