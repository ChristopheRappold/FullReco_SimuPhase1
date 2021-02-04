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

THyphiAttributes::THyphiAttributes(const FullRecoConfig& config, const DataSim& In)
  : Field(nullptr), Config(config), InputPar(In)
{
  _logger = spdlog::get("Console");
  _logger->info( "THyphiAttributes()" );

  
  G4_simu = false;
  G4_TimeResolution = false;
  G4_GeoResolution = false;
  Debug_DAF = false;
  DoNoMaterial = false;

  Task_CheckField = true;
  Task_FlatMCOutputML = false;
  Task_BayesFinder = false;
  Task_FinderCM = false;
  Task_CheckRZ = true;
  Task_KalmanDAF = true;

  RZ_ChangeMiniFiber = false;
  RZ_MDCProlate = true;
  RZ_MDCWire2 = false;
  RZ_MDCBiasCorr = true;

  KF_Kalman = false;
  KF_KalmanSqrt = true;
  KF_KalmanRef = false;
  KF_DAFRef = false;
  KF_DAF = false;

  KF_NbCentralCut = 1;
  KF_NbMiniFiberCut = 4;

  if(Config.IsAvailable("G4_simu"))
    G4_simu = true;
  if(Config.IsAvailable("G4_TimeReso"))
    G4_TimeResolution = true;
  if(Config.IsAvailable( "G4_GeoReso"))
    G4_GeoResolution = true;
  if(Config.IsAvailable( "back_tracking"))
    back_tracking = true;
  if(Config.IsAvailable( "beam_only"))
    beam_only = true;
  if(Config.IsAvailable( "Debug_DAF"))
    Debug_DAF = true;
  if(Config.IsAvailable( "NoMaterial"))
    DoNoMaterial = true;

  if(Config.IsAvailable("Task_CheckField"))
    Task_CheckField = Config.Get<bool>("Task_CheckField");
  if(Config.IsAvailable("Task_FlatMCOutputML"))
    Task_FlatMCOutputML = Config.Get<bool>("Task_FlatMCOutputML");;
  if(Config.IsAvailable("Task_BayesFinder"))
    Task_BayesFinder = Config.Get<bool>("Task_BayesFinder");
  if(Config.IsAvailable("Task_FinderCM"))
    Task_FinderCM = Config.Get<bool>("Task_FinderCM");
  if(Config.IsAvailable("Task_CheckRZ"))
    Task_CheckRZ = Config.Get<bool>("Task_CheckRZ");
  if(Config.IsAvailable("Task_KalmanDAF"))
    Task_KalmanDAF = Config.Get<bool>("Task_KalmanDAF");

  if(Config.IsAvailable("RZ_ChangeMiniFiber"))
    RZ_ChangeMiniFiber = Config.Get<bool>("RZ_ChangeMiniFiber");
  if(Config.IsAvailable("RZ_MDCProlate"))
    RZ_MDCProlate = Config.Get<bool>("RZ_MDCProlate");
  if(Config.IsAvailable("RZ_MDCWire2"))
    RZ_MDCWire2 = Config.Get<bool>("RZ_MDCWire2");
  if(Config.IsAvailable("RZ_MDCBiasCorr"))
    RZ_MDCBiasCorr = Config.Get<bool>("RZ_MDCBiasCorr");

  _logger->info("RZ Settting: Prolate? {} Wire2? {} BiasCorr? {} ChangeMiniF? {}",RZ_MDCProlate, RZ_MDCWire2, RZ_MDCBiasCorr, RZ_ChangeMiniFiber);

  if(Config.IsAvailable("KF_Kalman"))
    KF_Kalman = Config.Get<bool>("KF_Kalman");
  if(Config.IsAvailable("KF_KalmanRef"))
    KF_KalmanRef = Config.Get<bool>("KF_KalmanRef");
  if(Config.IsAvailable("KF_KalmanSqrt"))
    KF_KalmanSqrt = Config.Get<bool>("KF_KalmanSqrt");
  if(Config.IsAvailable("KF_DAF"))
    KF_DAF = Config.Get<bool>("KF_DAF");
  if(Config.IsAvailable("KF_DAFRef"))
    KF_DAFRef = Config.Get<bool>("KF_DAFRef");

  _logger->info("KF Setting: Kalman? {} / KalmanRef? {} / KalmanSqrt? {} / DAF? {} / DAFRef? {}",KF_Kalman, KF_KalmanRef, KF_KalmanSqrt, KF_DAF, KF_DAFRef);

  if(Config.IsAvailable("KF_NbCentralCut"))
    KF_NbCentralCut= Config.Get<int>("KF_NbCentralCut");
  if(Config.IsAvailable("KF_NbMiniFiberCut"))
    KF_NbMiniFiberCut= Config.Get<int>("KF_NbMiniFiberCut");

  _logger->info("KF_RejectionCut : Central < : {} MiniFiber < : {}",KF_NbCentralCut, KF_NbMiniFiberCut);

  std::string temp_name_out = config.Get<std::string>("Output_Namefile");
  std::string temp_file_base_name = temp_name_out.substr(0,temp_name_out.find_last_of('.'));

  std::string MLSuffix = Config.IsAvailable("FlatML_Suffix") ? Config.Get<std::string>("FlatML_Suffix") : "ML_MCOutput";

  temp_file_base_name += MLSuffix;
  temp_file_base_name += ".root";
  FlatML_namefile = temp_file_base_name;

  Nb_CPU = Config.Get<int>("Nb_CPU");
  Nb_Fraction = Config.Get<int>("Nb_Fraction");
  NEvent = Config.Get<int>("Nb_Event");
  if(Config.IsAvailable("MultiThreading"))
    {
      RunType = MultiTh;
      IsMain = true;
      NQueue   = Config.IsAvailable("Queue_Size")     ? Config.Get<int>("Queue_Size")     : 10;
      NBuilder = Config.IsAvailable("Nb_DataBuilder") ? Config.Get<int>("Nb_DataBuilder") : 1;
      NKalman  = Config.IsAvailable("Nb_TrackFitter") ? Config.Get<int>("Nb_TrackFitter") : 1;
      NMerger  = Config.IsAvailable("Nb_DataMerger")  ? Config.Get<int>("Nb_DataMerger")  : 1;
    }
  else if(Config.IsAvailable("ZeroMQ"))
    {
      RunType = ZeroMQ;
      NQueue   = -1;
      IsMain = Config.IsAvailable("ZMQ_Main") ? Config.Get<bool>("ZMQ_Main") : true;
      NBuilder = Config.IsAvailable("Nb_DataBuilder") ? Config.Get<int>("Nb_DataBuilder") : 1;
      NKalman  = Config.IsAvailable("Nb_TrackFitter") ? Config.Get<int>("Nb_TrackFitter") : 1;
      NMerger  = Config.IsAvailable("Nb_DataMerger")  ? Config.Get<int>("Nb_DataMerger")  : 1;

      addr_initEvent = Config.IsAvailable("Addr_InitEvent") ? Config.Get<std::string>("Addr_InitEvent") : "inproc://Nevent";
      
      addr_frontBuilder = Config.IsAvailable("Addr_InputBuilder") ? Config.Get<std::string>("Addr_InputBuilder") : "inproc://Q0";
      addr_backFitter = Config.IsAvailable("Addr_OutputBuilder") ? Config.Get<std::string>("Addr_OutputBuilder") : "inproc://Q0out";
      
      addr_frontFitter = Config.IsAvailable("Addr_InputFitter") ? Config.Get<std::string>("Addr_InputFitter") : "inproc://Q1";
      addr_backMerger = Config.IsAvailable("Addr_OutputFitter") ? Config.Get<std::string>("Addr_OutputFitter") : "inproc://Q1out";
      
      addr_frontMerger = Config.IsAvailable("Addr_InputMerger") ? Config.Get<std::string>("Addr_InputMerger") : "inproc://Q2";
      addr_backEnd = Config.IsAvailable("Addr_OutputMerger") ? Config.Get<std::string>("Addr_OutputMerger") : "inproc://Q2out";
      
      addr_control = Config.IsAvailable("Addr_EndControl") ? Config.Get<std::string>("Addr_EndControl") : "inproc://controlQ";

      addr_monitor = Config.IsAvailable("Addr_Monitor") ? Config.Get<std::string>("Addr_Monitor") : "tcp://127.0.0.1:9876";
    }
  else
    {
      RunType = SingleTh;
      NQueue   = -1;
      NBuilder = 1;
      NKalman  = 1;
      NMerger  = 1;
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


  gGeoManager->CloseGeometry();
  if(Config.IsAvailable("MultiThreading"))
    gGeoManager->SetMaxThreads(NKalman);
  
  _logger->info( " done ");
}

int THyphiAttributes::Init_Para()
{
  //Nb_CPU = 1;
  //Nb_Fraction = 1;

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
