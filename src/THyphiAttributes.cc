#include "THyphiAttributes.h"

#include <cstdlib>
#include <fstream>
#include <functional>
//#include <iostream>

#include "FairBase/WasaSolenoidFieldMap.h"
#include "TDatime.h"
#include "TMath.h"
#include "TTimeStamp.h"
//#define FIELD_DEBUG
#include "FieldManager.h"
#include "GFHypFieldMap_new.h"
#include "GFWasaMap.h"
#include "MaterialEffects.h"
#include "TGeoManager.h"
#include "TGeoMaterialInterface.h"
#include "spdlog/spdlog.h"

#include "boost/algorithm/string.hpp"

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

THyphiAttributes::THyphiAttributes(const FullRecoConfig& config, const DataSimExp& In)
    : Field(nullptr), Config(config), InputPar(In)
{
  _logger = spdlog::get("Console");
  _logger->info("THyphiAttributes()");

  NameIn  = config.Get<std::string>("Input_Namefile");
  NameOut = config.Get<std::string>("Output_Namefile");

  TTimeStamp timing;
  int timeA = timing.GetTime();
  int dateA = timing.GetDate();
  TDatime datetime;
  datetime.Set(dateA, timeA);

  DateOfRun = datetime.AsSQLString();

  std::hash<std::string> HashStr;
  std::string tempStr(NameIn);
  tempStr += NameOut;
  tempStr += DateOfRun;
  Hash = std::to_string(HashStr(tempStr));

  _logger->info("Id Hash : {} / dateOfRun : {}", Hash, DateOfRun);

  G4_simu           = false;
  G4_TimeResolution = false;
  G4_GeoResolution  = false;
  Debug_DAF         = false;
  DoNoMaterial      = false;

  RZ_ChangeMiniFiber = false;
  RZ_MDCProlate      = true;
  RZ_MDCWire2        = false;
  RZ_MDCBiasCorr     = true;

  KF_Kalman     = false;
  KF_KalmanSqrt = true;
  KF_KalmanRef  = false;
  KF_DAFRef     = false;
  KF_DAF        = false;

  KF_NbCentralCut   = 1;
  KF_NbMiniFiberCut = 4;

  RF_OutputEvents = false;

  TaskConfig.Init(Config);

  if(Config.IsAvailable("G4_simu"))
    G4_simu = true;
  if(Config.IsAvailable("G4_TimeReso"))
    G4_TimeResolution = true;
  if(Config.IsAvailable("G4_GeoReso"))
    G4_GeoResolution = true;
  if(Config.IsAvailable("back_tracking"))
    back_tracking = true;
  if(Config.IsAvailable("beam_only"))
    beam_only = true;
  if(Config.IsAvailable("Debug_DAF"))
    Debug_DAF = true;
  if(Config.IsAvailable("NoMaterial"))
    DoNoMaterial = true;
  if(Config.IsAvailable("RZ_ChangeMiniFiber"))
    RZ_ChangeMiniFiber = Config.Get<bool>("RZ_ChangeMiniFiber");
  if(Config.IsAvailable("RZ_MDCProlate"))
    RZ_MDCProlate = Config.Get<bool>("RZ_MDCProlate");
  if(Config.IsAvailable("RZ_MDCWire2"))
    RZ_MDCWire2 = Config.Get<bool>("RZ_MDCWire2");
  if(Config.IsAvailable("RZ_MDCBiasCorr"))
    RZ_MDCBiasCorr = Config.Get<bool>("RZ_MDCBiasCorr");

  _logger->info("RZ Settting: Prolate? {} Wire2? {} BiasCorr? {} ChangeMiniF? {}", RZ_MDCProlate, RZ_MDCWire2,
                RZ_MDCBiasCorr, RZ_ChangeMiniFiber);

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

  _logger->info("KF Setting: Kalman? {} / KalmanRef? {} / KalmanSqrt? {} / DAF? {} / DAFRef? {}", KF_Kalman,
                KF_KalmanRef, KF_KalmanSqrt, KF_DAF, KF_DAFRef);

  if(Config.IsAvailable("KF_NbCentralCut"))
    KF_NbCentralCut = Config.Get<int>("KF_NbCentralCut");
  if(Config.IsAvailable("KF_NbMiniFiberCut"))
    KF_NbMiniFiberCut = Config.Get<int>("KF_NbMiniFiberCut");

  _logger->info("KF_RejectionCut : Central < : {} MiniFiber < : {}", KF_NbCentralCut, KF_NbMiniFiberCut);

  StudyCase = Config.IsAvailable("StudyCase") ? Config.Get<std::string>("StudyCase") : "None";
  // std::string temp_name_out = config.Get<std::string>("Output_Namefile");
  std::string temp_file_base_name = NameOut.substr(0, NameOut.find_last_of('.'));

  std::string MLSuffix = Config.IsAvailable("FlatML_Suffix") ? Config.Get<std::string>("FlatML_Suffix") : "ML_MCOutput";
  temp_file_base_name += MLSuffix;
  temp_file_base_name += ".root";
  FlatML_namefile = temp_file_base_name;

  if(Config.IsAvailable("RF_OutputEvents"))
    RF_OutputEvents = Config.Get<bool>("RF_OutputEvents");

  DataML_Out = Config.IsAvailable("DataML_Out") ? Config.Get<std::string>("DataML_Out") : "NoneInConfig";

  Nb_CPU      = Config.Get<int>("Nb_CPU");
  Nb_Fraction = Config.Get<int>("Nb_Fraction");
  NEvent      = Config.Get<int>("Nb_Event");

  MTsetting.Init(Config);

  if(TaskConfig.Task_ReStart)
    {
      int res = Reload_Para();
      if(res != 0)
	_logger->info(" E> Reload config did not happen correctly ! ");
      else
	_logger->info(" *** > Reload config done ! ");
    }
  else
    Init_Para();

  _logger->info(" *** > Loading Fieldmap ");

if(Config.IsAvailable("CalibFile_MDCmap"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_MDCmap");
    map_ParamFiles.insert({"mdc_map_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"mdc_map_file","./calib/mdc_mapping/MDC_channelmap.csv"});

if(Config.IsAvailable("CalibFile_MDCphys"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_MDCphys");
    map_ParamFiles.insert({"mdc_phys_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"mdc_phys_file","./calib/mdc_mapping/MDC_PhysicalMap.csv"});

if(Config.IsAvailable("CalibFile_MDCpar"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_MDCpar");
    map_ParamFiles.insert({"mdc_par_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"mdc_par_file","./calib/mdc_driftparam/MDC_DriftParam.txt"});

if(Config.IsAvailable("CalibFile_Fiberoffset"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_Fiberoffset");
    map_ParamFiles.insert({"fiber_offset_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"fiber_offset_file","./calib/fiber_offset/fiber_offset.csv"});

if(Config.IsAvailable("CalibFile_PSBtime"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_PSBtime");
    map_ParamFiles.insert({"psb_time_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"psb_time_file","calib/psb_param/psb_time.csv"});

if(Config.IsAvailable("CalibFile_T0time"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_T0time");
    map_ParamFiles.insert({"t0_time_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"t0_time_file","./calib/t0_param/t0_time.csv"});



  if(Wasa_FieldMap)
    {
      double signDir = Wasa_Side ? 1.0 : -1.0;
      Field = new WasaSolenoidFieldMap("WasaFieldMap", "WasaFieldMap", Wasa_FieldMapName, Field_Strength, signDir);
    }
  else
    Field = new FrsSolenoidHypField();

  bool isWasa = false;
  for(const auto& nameGeo : name_GeoVolumes)
    if(nameGeo == "WASA")
      {
        isWasa = true;
        _logger->warn("!> Wasa geometry found !");
        break;
      }

  if(isWasa == false)
    dynamic_cast<FrsSolenoidHypField*>(Field)->SetPositionFromGeoManager("CDS_logR_0");
  else
    {
      if(Wasa_FieldMap)
        dynamic_cast<WasaSolenoidFieldMap*>(Field)->SetPositionFromGeoManager("MFLD_1");
      else
        dynamic_cast<FrsSolenoidHypField*>(Field)->SetPositionFromGeoManager("INNER_1");
    }

  if(Wasa_FieldMap)
    {
      dynamic_cast<WasaSolenoidFieldMap*>(Field)->Init();
      dynamic_cast<WasaSolenoidFieldMap*>(Field)->Print();
    }
  else
    {
      dynamic_cast<FrsSolenoidHypField*>(Field)->SetField(0., 0., Field_Strength);
      dynamic_cast<FrsSolenoidHypField*>(Field)->Print();
    }
  // std::cout<< " IsInside : (0,0,70)"<<  dynamic_cast<FrsSolenoidHypField*>(Field)->IsInside(0.,0.,70.)<< "
  // (0,0,125.50) "<<  dynamic_cast<FrsSolenoidHypField*>(Field)->IsInside(0.,0.,125.5)<<"\n"; std::cout<< " IsInside :
  // (0,0,40)"<<  dynamic_cast<FrsSolenoidHypField*>(Field)->IsInside(0.,0.,40.)<< " (0,0,155.50) "<<
  // dynamic_cast<FrsSolenoidHypField*>(Field)->IsInside(0.,0.,155.5)<<"\n";

  // bool inTelsa = false;
  // GFFieldManager::getInstance()->init(new GFHypFieldMap_new(true,-0.750533/0.8016*facFRS,inTelsa,false,1.,Field));
  // double facFRS = Field_Strength;
  // genfit::FieldManager::getInstance()->init(new genfit::GFHypFieldMap_new(true, 0.750533 * facFRS, inTelsa,
  // false, 1., false, Field));
  genfit::FieldManager::getInstance()->init(new genfit::GFWasaMap(Field));

  TGeoMedium* vac = gGeoManager->GetMedium(" VAC");
  assert(vac != nullptr);
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
      _logger->warn(" ** > Use No material !");
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
    gGeoManager->SetMaxThreads(MTsetting.NKalman);

  _logger->info(" done ");
}

int THyphiAttributes::Init_Para()
{
  // Nb_CPU = 1;
  // Nb_Fraction = 1;

  auto posTargetX  = InputPar.simParameters->find("Target_PosX");
  Target_PositionX = posTargetX->second;
  auto posTargetY  = InputPar.simParameters->find("Target_PosY");
  Target_PositionY = posTargetY->second;
  auto posTargetZ  = InputPar.simParameters->find("Target_PosZ");
  Target_PositionZ = posTargetZ->second;
  auto sizeTarget  = InputPar.simParameters->find("Target_Size");
  Target_Size      = sizeTarget->second;
  auto fieldV      = InputPar.simParameters->find("Field_CDS_Bz");
  Field_Strength   = fieldV->second;
  auto WasaS       = InputPar.simParameters->find("Wasa_Side");
  Wasa_Side        = WasaS->second;

  auto WasaMap  = InputPar.simParameters->find("Field_CDS_FieldMap");
  Wasa_FieldMap = WasaMap != InputPar.simParameters->end() ? true : false;
  if(Wasa_FieldMap)
    Wasa_FieldMapName = "./field/MagField_default.dat";

  TObjArray* L_vol = gGeoManager->GetListOfVolumes();
  int n_volume     = L_vol->GetEntries();
  for(int n_v = 0; n_v < n_volume; ++n_v)
    {
      std::string name_vol_temp(L_vol->At(n_v)->GetName());
      name_GeoVolumes.push_back(name_vol_temp);
    }
  name_GeoVolumes.push_back("Total");

  return 0;
}

int THyphiAttributes::Reload_Para()
{
  std::string Hash(InputPar.previousMeta->Hash);

  auto storage = InitStorage();
  storage.sync_schema();

  auto previousConfSel = storage.get_all<RunAttrGeneralDef>(sqlite_orm::where(sqlite_orm::is_equal(&RunAttrGeneralDef::Hash,Hash)));

  if(previousConfSel.size()!=1)
    {
      _logger->info("Reload from database : could not found the proper previous configuration : {} : nb Entries {}",Hash,previousConfSel.size());
      return -1;
    }

  auto previousConf = previousConfSel.front();

  Target_PositionX = previousConf.Target_PositionX;
  Target_PositionY = previousConf.Target_PositionY;
  Target_PositionZ = previousConf.Target_PositionZ;
  Target_Size      = previousConf.Target_Size;
  Field_Strength   = previousConf.Field_Strength;
  Wasa_Side        = previousConf.Wasa_Side;

  Wasa_FieldMap = previousConf.Wasa_FieldMap;
  Wasa_FieldMapName = previousConf.Wasa_FieldMapName;

  std::string tempStr;
  for(size_t i=0; i<previousConf.name_GeoVolumes.size();++i)
    {
      if(previousConf.name_GeoVolumes[i] != ';')
	tempStr += previousConf.name_GeoVolumes[i];
      else
	{
	  name_GeoVolumes.push_back(tempStr);
	  tempStr.clear();
	}
    }

  auto previousFullConfSel = storage.get_all<RunFullConfigDef>(sqlite_orm::where(sqlite_orm::is_equal(&RunFullConfigDef::Hash,Hash)));

  auto previousFullConf = previousFullConfSel.front();

  int diffConf = Config.Reload(previousFullConf.ConfigJson);
  if(diffConf !=0)
    return -2;

  return 0;
}


void MT::Init(const FullRecoConfig& Config)
{
  if(Config.IsAvailable("MultiThreading"))
    {
      RunType  = MultiTh;
      IsMain   = true;
      NQueue   = Config.IsAvailable("Queue_Size") ? Config.Get<int>("Queue_Size") : 10;
      NBuilder = Config.IsAvailable("Nb_DataBuilder") ? Config.Get<int>("Nb_DataBuilder") : 1;
      NKalman  = Config.IsAvailable("Nb_TrackFitter") ? Config.Get<int>("Nb_TrackFitter") : 1;
      NMerger  = Config.IsAvailable("Nb_DataMerger") ? Config.Get<int>("Nb_DataMerger") : 1;
    }
  else if(Config.IsAvailable("ZeroMQ"))
    {
      RunType  = ZeroMQ;
      NQueue   = -1;
      IsMain   = Config.IsAvailable("ZMQ_Main") ? Config.Get<bool>("ZMQ_Main") : true;
      NBuilder = Config.IsAvailable("Nb_DataBuilder") ? Config.Get<int>("Nb_DataBuilder") : 1;
      NKalman  = Config.IsAvailable("Nb_TrackFitter") ? Config.Get<int>("Nb_TrackFitter") : 1;
      NMerger  = Config.IsAvailable("Nb_DataMerger") ? Config.Get<int>("Nb_DataMerger") : 1;

      addr_initEvent =
          Config.IsAvailable("Addr_InitEvent") ? Config.Get<std::string>("Addr_InitEvent") : "inproc://Nevent";

      addr_frontBuilder =
          Config.IsAvailable("Addr_InputBuilder") ? Config.Get<std::string>("Addr_InputBuilder") : "inproc://Q0";
      addr_backFitter =
          Config.IsAvailable("Addr_OutputBuilder") ? Config.Get<std::string>("Addr_OutputBuilder") : "inproc://Q0out";

      addr_frontFitter =
          Config.IsAvailable("Addr_InputFitter") ? Config.Get<std::string>("Addr_InputFitter") : "inproc://Q1";
      addr_backMerger =
          Config.IsAvailable("Addr_OutputFitter") ? Config.Get<std::string>("Addr_OutputFitter") : "inproc://Q1out";

      addr_frontMerger =
          Config.IsAvailable("Addr_InputMerger") ? Config.Get<std::string>("Addr_InputMerger") : "inproc://Q2";
      addr_backEnd =
          Config.IsAvailable("Addr_OutputMerger") ? Config.Get<std::string>("Addr_OutputMerger") : "inproc://Q2out";

      addr_control =
          Config.IsAvailable("Addr_EndControl") ? Config.Get<std::string>("Addr_EndControl") : "inproc://controlQ";

      addr_monitor =
          Config.IsAvailable("Addr_Monitor") ? Config.Get<std::string>("Addr_Monitor") : "tcp://127.0.0.1:9876";
    }
  else
    {
      RunType  = SingleTh;
      NQueue   = -1;
      NBuilder = 1;
      NKalman  = 1;
      NMerger  = 1;
    }
}

void Task::Init(const FullRecoConfig& Config)
{
  if(Config.IsAvailable("Task_CheckField"))
    Task_CheckField = Config.Get<bool>("Task_CheckField");
  if(Config.IsAvailable("Task_PrimaryVtx"))
    Task_PrimaryVtx = Config.Get<bool>("Task_PrimaryVtx");
  if(Config.IsAvailable("Task_FlatMCOutputML"))
    Task_FlatMCOutputML = Config.Get<bool>("Task_FlatMCOutputML");
  if(Config.IsAvailable("Task_CheckFiberTrack"))
    Task_CheckFiberTrack = Config.Get<bool>("Task_CheckFiberTrack");
  if(Config.IsAvailable("Task_BayesFinder"))
    Task_BayesFinder = Config.Get<bool>("Task_BayesFinder");
  if(Config.IsAvailable("Task_RiemannFinder"))
    Task_RiemannFinder = Config.Get<bool>("Task_RiemannFinder");
  if(Config.IsAvailable("Task_FindingPerf"))
    Task_FindingPerf = Config.Get<bool>("Task_FindingPerf");
  if(Config.IsAvailable("Task_FinderCM"))
    Task_FinderCM = Config.Get<bool>("Task_FinderCM");
  if(Config.IsAvailable("Task_CheckRZ"))
    Task_CheckRZ = Config.Get<bool>("Task_CheckRZ");
  if(Config.IsAvailable("Task_KalmanDAF"))
    Task_KalmanDAF = Config.Get<bool>("Task_KalmanDAF");
  if(Config.IsAvailable("Task_DecayVtx"))
    Task_DecayVtx = Config.Get<bool>("Task_DecayVtx");
  if(Config.IsAvailable("Task_ReStart"))
    Task_ReStart = Config.Get<bool>("Task_ReStart");


  if(Config.IsAvailable("Task_Pipeline"))
    {
      Task_Order.clear();

      std::string in(Config.Get<std::string>("Task_Pipeline"));
      std::vector<std::string> pipeline;
      boost::split(pipeline, in, boost::is_any_of("->"), boost::token_compress_on);

      for(auto& s : pipeline)
	boost::trim(s);

      for(const auto& s : pipeline)
	{
	  if(s == "Task_CheckField")
	    Task_Order.push_back(TASKCHECKFIELD);
	  if(s == "Task_PrimaryVtx")
	    Task_Order.push_back(TASKPRIMARYVTX);
	  if(s == "Task_FlatMCOutputML")
	    Task_Order.push_back(TASKFLATMCOUTPUTML);
	  if(s == "Task_CheckFiberTrack")
	    Task_Order.push_back(TASKCHECKFIBERTRACK);
  	  if(s == "Task_BayesFinder")
	    Task_Order.push_back(TASKBAYESFINDER);
	  if(s == "Task_RiemannFinder")
	    Task_Order.push_back(TASKRIEMANNFINDER);
	  if(s == "Task_FindingPerf")
	    Task_Order.push_back(TASKFINDINGPERF);
	  if(s == "Task_FinderCM")
	    Task_Order.push_back(TASKFINDERCM);
	  if(s == "Task_CheckRZ")
	    Task_Order.push_back(TASKCHECKRZ);
	  if(s == "Task_KalmanDAF")
	    Task_Order.push_back(TASKKALMANDAF);
	  if(s == "Task_DecayVtx")
	    Task_Order.push_back(TASKDECAYVTX);
	  if(s == "Task_ReStart")
	    Task_Order.push_back(TASKRESTART);
	}
    }
}

void THyphiAttributes::SetOut(AttrOut& out) const
{
  out.RunDone.Hash      = Hash;
  out.RunDone.NameIn    = NameIn;
  out.RunDone.NameOut   = NameOut;
  out.RunDone.DateOfRun = DateOfRun;
  out.RunDone.Nb_CPU               = Nb_CPU;
  out.RunDone.Nb_Fraction          = Nb_Fraction;
  out.RunDone.NEvent               = NEvent;

  out.RunMulti.Hash = Hash;
  out.RunMulti.MT_RunType           = MTsetting.RunType;
  out.RunMulti.MT_IsMain            = MTsetting.IsMain;
  out.RunMulti.MT_NQueue            = MTsetting.NQueue;
  out.RunMulti.MT_NBuilder          = MTsetting.NBuilder;
  out.RunMulti.MT_NKalman           = MTsetting.NKalman;
  out.RunMulti.MT_NMerger           = MTsetting.NMerger;
  out.RunMulti.MT_addr_initEvent    = MTsetting.addr_initEvent;
  out.RunMulti.MT_addr_frontBuilder = MTsetting.addr_frontBuilder;
  out.RunMulti.MT_addr_backFitter   = MTsetting.addr_backFitter;
  out.RunMulti.MT_addr_frontFitter  = MTsetting.addr_frontFitter;
  out.RunMulti.MT_addr_backMerger   = MTsetting.addr_backMerger;
  out.RunMulti.MT_addr_frontMerger  = MTsetting.addr_frontMerger;
  out.RunMulti.MT_addr_backEnd      = MTsetting.addr_backEnd;
  out.RunMulti.MT_addr_control      = MTsetting.addr_control;
  out.RunMulti.MT_addr_monitor      = MTsetting.addr_monitor;

  out.RunTask.Hash = Hash;
  out.RunTask.Task_CheckField      = TaskConfig.Task_CheckField;
  out.RunTask.Task_PrimaryVtx      = TaskConfig.Task_PrimaryVtx;
  out.RunTask.Task_FlatMCOutputML  = TaskConfig.Task_FlatMCOutputML;
  out.RunTask.Task_BayesFinder     = TaskConfig.Task_BayesFinder;
  out.RunTask.Task_RiemannFinder   = TaskConfig.Task_RiemannFinder;
  out.RunTask.Task_FinderCM        = TaskConfig.Task_FinderCM;
  out.RunTask.Task_FindingPerf     = TaskConfig.Task_FindingPerf;
  out.RunTask.Task_CheckRZ         = TaskConfig.Task_CheckRZ;
  out.RunTask.Task_KalmanDAF       = TaskConfig.Task_KalmanDAF;
  out.RunTask.Task_DecayVtx        = TaskConfig.Task_DecayVtx;
  out.RunTask.Task_ReStart         = TaskConfig.Task_ReStart;

  out.RunAttrGeneral.Hash = Hash;
  out.RunAttrGeneral.G4_simu              = G4_simu;
  out.RunAttrGeneral.G4_TimeResolution    = G4_TimeResolution;
  out.RunAttrGeneral.G4_GeoResolution     = G4_GeoResolution;
  out.RunAttrGeneral.back_tracking        = back_tracking;
  out.RunAttrGeneral.Target_PositionX     = Target_PositionX;
  out.RunAttrGeneral.Target_PositionY     = Target_PositionY;
  out.RunAttrGeneral.Target_PositionZ     = Target_PositionZ;
  out.RunAttrGeneral.Target_Size          = Target_Size;
  out.RunAttrGeneral.Field_Strength       = Field_Strength;
  out.RunAttrGeneral.Wasa_Side            = Wasa_Side;
  out.RunAttrGeneral.Wasa_FieldMap        = Wasa_FieldMap;
  out.RunAttrGeneral.Wasa_FieldMapName    = Wasa_FieldMapName;

  for(auto nameG : name_GeoVolumes)
    {
      for(size_t i = 0; i < nameG.length(); ++i)
        out.RunAttrGeneral.name_GeoVolumes.push_back(nameG[i]);
      out.RunAttrGeneral.name_GeoVolumes.push_back(';');
    }

  out.RunAttrGeneral.beam_only          = beam_only;

  out.RunTaskAttr.Hash = Hash;
  out.RunTaskAttr.Debug_DAF          = Debug_DAF;
  out.RunTaskAttr.DoNoMaterial       = DoNoMaterial;
  out.RunTaskAttr.RZ_ChangeMiniFiber = RZ_ChangeMiniFiber;
  out.RunTaskAttr.RZ_MDCProlate      = RZ_MDCProlate;
  out.RunTaskAttr.RZ_MDCWire2        = RZ_MDCWire2;
  out.RunTaskAttr.RZ_MDCBiasCorr     = RZ_MDCBiasCorr;
  out.RunTaskAttr.KF_Kalman          = KF_Kalman;
  out.RunTaskAttr.KF_KalmanSqrt      = KF_KalmanSqrt;
  out.RunTaskAttr.KF_KalmanRef       = KF_KalmanRef;
  out.RunTaskAttr.KF_DAFRef          = KF_DAFRef;
  out.RunTaskAttr.KF_DAF             = KF_DAF;
  out.RunTaskAttr.KF_NbCentralCut    = KF_NbCentralCut;
  out.RunTaskAttr.KF_NbMiniFiberCut  = KF_NbMiniFiberCut;
  out.RunTaskAttr.FlatML_namefile    = FlatML_namefile;
  out.RunTaskAttr.DataML_Out         = DataML_Out;
  out.RunTaskAttr.RF_OutputEvents    = RF_OutputEvents;

  out.RunFullConfig.Hash = Hash;
  std::string tempCj;
  Config.Save(tempCj);
  for(size_t i = 0; i < tempCj.length(); ++i)
    out.RunFullConfig.ConfigJson.push_back(tempCj[i]);
}

void THyphiAttributes::SaveToDatabase() const
{
  auto storage = InitStorage();
  storage.sync_schema(true);

  AttrOut outA;
  SetOut(outA);
  _logger->debug("ToDatabase : hash {} ", outA.RunDone.Hash);
  // storage.begin_transaction();
  storage.insert(outA.RunDone);
  storage.insert(outA.RunMulti);
  storage.insert(outA.RunTask);
  storage.insert(outA.RunAttrGeneral);
  storage.insert(outA.RunTaskAttr);
  storage.insert(outA.RunFullConfig);
  // storage.commit();
}
