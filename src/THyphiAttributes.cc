#include "THyphiAttributes.h"
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "TMath.h"
//#define FIELD_DEBUG
#include "FieldManager.h"
#include "GFHypFieldMap_new.h"
#include "GFWasaMap.h"
#include "MaterialEffects.h"
#include "TGeoMaterialInterface.h"

#include "TGeoManager.h"

//#define OLDGENFIT For Kalman_RK
//#define OLD_FIELDMAP
THyphiAttributes::THyphiAttributes() : Field(nullptr), InputPar{nullptr, nullptr}
{

  G4_simu = false;
  G4_TimeResolution = false;
  G4_GeoResolution = false;
  back_tracking = false;
  beam_only = false;
  Debug_DAF = false;
  DoNoMaterial = false;

  Init_Para();
}

THyphiAttributes::THyphiAttributes(const std::list<std::string>& type, const std::list<std::string>& option, double FieldScalingFactor,
                                   const DataSim& In)
    : Field(nullptr), InputPar(In)
{
  std::cout << "THyphiAttributes()" << std::endl;

  G4_simu = false;
  G4_TimeResolution = false;
  G4_GeoResolution = false;
  Debug_DAF = false;
  DoNoMaterial = false;

  for(const auto& opt : option)
    {
      if(opt == "G4_simu")
        G4_simu = true;
      else if(opt == "G4_TimeReso")
        G4_TimeResolution = true;
      else if(opt == "G4_GeoReso")
        G4_GeoResolution = true;
      else if(opt == "back_tracking")
        back_tracking = true;
      else if(opt == "beam_only")
        beam_only = true;
      else if(opt == "Debug_DAF")
        Debug_DAF = true;
      else if(opt == "NoMaterial")
        DoNoMaterial = true;
      else
        {
          size_t found;
          found = opt.find("CPU");
          if(found != std::string::npos)
            {
              std::string substring(opt.substr(found + 3));
              Nb_CPU = std::atoi(substring.c_str());
            }
          found = opt.find("Fraction");
          if(found != std::string::npos)
            {
              std::string substring(opt.substr(found + 8));
              Nb_Fraction = std::atoi(substring.c_str());
            }
	  found = opt.find("Event");
	  if(found != std::string::npos)
	    {
	      std::string substring(opt.substr(found + 5));
              NEvent = std::atoi(substring.c_str());
	    }
        }
    }

  Init_Para();

  std::cout << " *** > Loading Fieldmap ";

  Field = new FrsSolenoidHypField();

  bool isWasa = false;
  for(auto nameGeo : name_GeoVolumes)
    if(nameGeo == "WASA")
      {
        isWasa = true;
        std::cout << "!> Wasa geometry found !\n";
        break;
      }

  if(isWasa == false)
    dynamic_cast<FrsSolenoidHypField*>(Field)->SetPositionFromGeoManager("CDS_logR_0");
  else
    dynamic_cast<FrsSolenoidHypField*>(Field)->SetPositionFromGeoManager("INNER_1");

  dynamic_cast<FrsSolenoidHypField*>(Field)->SetField(0.,0.,Field_Strength);
  //bool inTelsa = false;
  // GFFieldManager::getInstance()->init(new GFHypFieldMap_new(true,-0.750533/0.8016*facFRS,inTelsa,false,1.,Field));
  //double facFRS = Field_Strength;
  //genfit::FieldManager::getInstance()->init(new genfit::GFHypFieldMap_new(true, 0.750533 * facFRS, inTelsa, false, 1., false, Field));
  genfit::FieldManager::getInstance()->init(new genfit::GFWasaMap(Field));

  genfit::FieldManager::getInstance()->useCache(true, 8);
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  if(DoNoMaterial)
    {
      genfit::MaterialEffects::getInstance()->setNoEffects();
      std::cout << " ** > Use No material !" << std::endl;
    }

  std::cout << " done " << std::endl;
}

int THyphiAttributes::Init_Para()
{
  Nb_CPU = 1;
  Nb_Fraction = 1;

  auto posTargetX = InputPar.simParameters->find("Target_PosX");
  Target_PositionX = posTargetX->second;
  auto posTargetY = InputPar.simParameters->find("Target_PosY");
  Target_PositionY = posTargetY->second;
  auto posTargetZ = InputPar.simParameters->find("Target_PosZ");
  Target_PositionZ = posTargetZ->second;
  auto sizeTarget = InputPar.simParameters->find("Target_Size");
  Target_Size = sizeTarget->second;
  auto fieldV = InputPar.simParameters->find("Field_CDS_Bz");
  Field_Strength = fieldV->second;

  TObjArray* L_vol = gGeoManager->GetListOfVolumes();
  int n_volume = L_vol->GetEntries();
  for(int n_v = 0; n_v < n_volume; ++n_v)
    {
      std::string name_vol_temp(L_vol->At(n_v)->GetName());
      name_GeoVolumes.push_back(name_vol_temp);
    }
  name_GeoVolumes.push_back("Total");

  return 0;
}

THyphiAttributes::~THyphiAttributes() {}
