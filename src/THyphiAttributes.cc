#include "THyphiAttributes.h"
#include <cstdlib>
#include <fstream>
//#include <iostream>

#include "TMath.h"
//#define FIELD_DEBUG
#include "FieldManager.h"
#include "GFHypFieldMap_new.h"
#include "GFWasaMap.h"
#include "MaterialEffects.h"
#include "TGeoMaterialInterface.h"

#include "TGeoManager.h"

#include "spdlog/spdlog.h"

//#define OLDGENFIT For Kalman_RK
//#define OLD_FIELDMAP
// THyphiAttributes::THyphiAttributes() : Field(nullptr), InputPar{nullptr, nullptr}
// {
//   _logger = spdlog::get("Console");
  
//   G4_simu = false;
//   G4_TimeResolution = false;
//   G4_GeoResolution = false;
//   back_tracking = false;
//   beam_only = false;
//   Debug_DAF = false;
//   DoNoMaterial = false;

//   Init_Para();
// }

THyphiAttributes::THyphiAttributes(const std::list<std::string>& type, const std::list<std::string>& option, double FieldScalingFactor,
                                   const DataSim& In)
    : Field(nullptr), InputPar(In)
{
  _logger = spdlog::get("Console");
  _logger->info( "THyphiAttributes()" );

  
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

  _logger->info( " *** > Loading Fieldmap ");

  Field = new FrsSolenoidHypField();

  bool isWasa = false;
  for(const auto& nameGeo : name_GeoVolumes)
    if(nameGeo == "WASA")
      {
        isWasa = true;
        _logger->warn( "!> Wasa geometry found !");
        break;
      }

  if(isWasa == false)
    dynamic_cast<FrsSolenoidHypField*>(Field)->SetPositionFromGeoManager("CDS_logR_0");
  else
    dynamic_cast<FrsSolenoidHypField*>(Field)->SetPositionFromGeoManager("INNER_1");

  dynamic_cast<FrsSolenoidHypField*>(Field)->SetField(0.,0.,Field_Strength);
  dynamic_cast<FrsSolenoidHypField*>(Field)->Print();
  // std::cout<< " IsInside : (0,0,70)"<<  dynamic_cast<FrsSolenoidHypField*>(Field)->IsInside(0.,0.,70.)<< " (0,0,125.50) "<<  dynamic_cast<FrsSolenoidHypField*>(Field)->IsInside(0.,0.,125.5)<<"\n";
  // std::cout<< " IsInside : (0,0,40)"<<  dynamic_cast<FrsSolenoidHypField*>(Field)->IsInside(0.,0.,40.)<< " (0,0,155.50) "<<  dynamic_cast<FrsSolenoidHypField*>(Field)->IsInside(0.,0.,155.5)<<"\n";

  //bool inTelsa = false;
  // GFFieldManager::getInstance()->init(new GFHypFieldMap_new(true,-0.750533/0.8016*facFRS,inTelsa,false,1.,Field));
  //double facFRS = Field_Strength;
  //genfit::FieldManager::getInstance()->init(new genfit::GFHypFieldMap_new(true, 0.750533 * facFRS, inTelsa, false, 1., false, Field));
  genfit::FieldManager::getInstance()->init(new genfit::GFWasaMap(Field));

  TGeoMedium* vac = gGeoManager->GetMedium(" VAC");
  assert(vac!=nullptr);
  gGeoManager->GetVolume("PSCE")->SetMedium(vac);
  gGeoManager->GetVolume("PSFE")->SetMedium(vac);
  gGeoManager->GetVolume("HypHI_RPC_l_log")->SetMedium(vac);
  gGeoManager->GetVolume("HypHI_RPC_h_log")->SetMedium(vac);
  gGeoManager->GetVolume("FMF2_log")->SetMedium(vac);
  //  for(auto name : name_GeoVolumes)
  //    {
  //      TGeoVolume* vol = gGeoManager->GetVolume(name.c_str());
  //      if(vol != nullptr)
  // 	{
  // 	  vol->SetMedium(vac);
  // 	}
  //   }
  
  
  genfit::FieldManager::getInstance()->useCache(true, 8);
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  if(DoNoMaterial)
    {
      genfit::MaterialEffects::getInstance()->setNoEffects();
      _logger->warn( " ** > Use No material !" );
    }

  // genfit::MaterialEffects::getInstance()->drawdEdx(-211);
  // genfit::MaterialEffects::getInstance()->drawdEdx(2212);
  // genfit::MaterialEffects::getInstance()->drawdEdx(10003);
  // double bx=0.,by=0.,bz=0.;
  // genfit::FieldManager::getInstance()->getFieldVal(0.,0.,70.,bx,by,bz);
  // std::cout<<"From genfit: (0,0,70) bz:"<<bz<<" ";
  // genfit::FieldManager::getInstance()->getFieldVal(0.,0.,125.5,bx,by,bz);
  // std::cout<<"(0,0,125.5) bz:"<<bz<<"\n";
  // genfit::FieldManager::getInstance()->getFieldVal(0.,0.,40.,bx,by,bz);
  // std::cout<<"From genfit: (0,0,40) bz:"<<bz<<" ";
  // genfit::FieldManager::getInstance()->getFieldVal(0.,0.,155.5,bx,by,bz);
  // std::cout<<"(0,0,155.5) bz:"<<bz<<"\n";

  
  _logger->info( " done ");
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
  auto WasaS = InputPar.simParameters->find("Wasa_Side");
  Wasa_Side = WasaS->second;

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
