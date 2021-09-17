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

}
