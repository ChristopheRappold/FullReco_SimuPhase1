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

  PV_RealXUVComb   = false;
  PV_RealPrimTrack = false;

  RZ_ChangeMiniFiber = false;
  RZ_MDCProlate      = true;
  RZ_MDCWire2        = false;
  RZ_MDCBiasCorr     = true;

  WASAFinder_perfect = false;
  WASAFinder_PSBHits = true;
  WASAFinder_PSFEHits= false;

  KF_Kalman     = false;
  KF_KalmanSqrt = true;
  KF_KalmanRef  = false;
  KF_DAFRef     = false;
  KF_DAF        = false;

  KF_NbCentralCut   = 1;
  KF_NbMiniFiberCut = 4;

  RF_OutputEvents = false;


  //Optics parameters
  if(Config.IsAvailable("Optics_name")) optics_name = Config.Get<std::string>("Optics_name");
  else                                  optics_name = "./calib/optics/optics_par.csv";

  std::ifstream ifs_optics ( optics_name );
  if(ifs_optics.is_open()){
    const std::string CommentSymbol("#");

    std::string temp_line;
    while(std::getline(ifs_optics,temp_line)){
      std::stringstream stream(temp_line);
      std::string testComment(stream.str());
      std::size_t it_comment = testComment.find(CommentSymbol);
      if(it_comment!=std::string::npos){
        //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
        continue;
      }
      char name[8];
      float par1;
      float par2;
      stream >> name >> par1  >> par2;

      optics_par[name].emplace_back(par1);
      optics_par[name].emplace_back(par2);
    }

    printf("Optics done : %s\n", optics_name.c_str());
  }
  else
  {
    printf(" ! fail to open  : %s\n", optics_name.c_str());
    exit(-1);
  }

  optics_s2z = 4187.5; //Ti: 4150 mm, take this value just for easy expression
  
  fiber_mft1_pos_z = 226.93 - 0.6 - 0.02;
  fiber_mft2_pos_z = 230.93 - 0.6 + 0.02;

  psb_pos_x = 0.5; psb_pos_y = 5.5; psb_pos_z = 2760.;
  psb_rot_z = -0.4;
  cut_psb_phi = 0.4;
  cut_psb_z   = 150;
  cut_phi_fm  = 0.3;

  flag_dup_trackhit = true;
  flag_dup_trackhit_mdc = false;
  flag_trackhit_inclusive = true;


  // PID Cuts
  if(Config.IsAvailable("Pion_Cut_name"))
    pi_cut_name = Config.Get<std::string>("Pion_Cut_name");
  else
    pi_cut_name = "./calib/pid_cuts/cut_pi.root";

  TFile *f_cut = new TFile(pi_cut_name.c_str());
  cut_pi = (TCutG*)f_cut->Get("cut_pi");

  printf("pion cut loaded : %s\n", pi_cut_name.c_str());

  // MWDC Par load
  if(Config.IsAvailable("MWDC_Par_FitorHist"))
    mwdc_par_fitorhist = Config.Get<int>("MWDC_Par_FitorHist");
  else
    mwdc_par_fitorhist = 0;


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

  if(Config.IsAvailable("PV_RealXUVComb"))
    PV_RealXUVComb = Config.Get<bool>("PV_RealXUVComb");
  if(Config.IsAvailable("PV_RealPrimTrack"))
    PV_RealPrimTrack = Config.Get<bool>("PV_RealPrimTrack");

  if(Config.IsAvailable("RZ_ChangeMiniFiber"))
    RZ_ChangeMiniFiber = Config.Get<bool>("RZ_ChangeMiniFiber");
  if(Config.IsAvailable("RZ_MDCProlate"))
    RZ_MDCProlate = Config.Get<bool>("RZ_MDCProlate");
  if(Config.IsAvailable("RZ_MDCWire2"))
    RZ_MDCWire2 = Config.Get<bool>("RZ_MDCWire2");
  if(Config.IsAvailable("RZ_MDCBiasCorr"))
    RZ_MDCBiasCorr = Config.Get<bool>("RZ_MDCBiasCorr");

  if(Config.IsAvailable("WASAFinder_perfect"))
    WASAFinder_perfect = Config.Get<bool>("WASAFinder_perfect");
  if(Config.IsAvailable("WASAFinder_PSBHits"))
    WASAFinder_PSBHits = Config.Get<bool>("WASAFinder_PSBHits");
  if(Config.IsAvailable("WASAFinder_PSFEHits"))
    WASAFinder_PSFEHits = Config.Get<bool>("WASAFinder_PSFEHits");
  
  _logger->info("RZ Setting: Prolate? {} / Wire2? {} / BiasCorr? {} / ChangeMiniF? {}", RZ_MDCProlate, RZ_MDCWire2,
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

if(Config.IsAvailable("CalibFile_MDCdrift"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_MDCdrift");
    map_ParamFiles.insert({"mdc_drift_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"mdc_drift_file","./calib/mdc_driftparam/MDC_DriftParam.txt"});

if(Config.IsAvailable("CalibFile_MDCparT0"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_MDCparT0");
    map_ParamFiles.insert({"mdc_parT0_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"mdc_parT0_file","./calib/mdc_driftparam/MDC_T0Param.csv"});

if(Config.IsAvailable("CalibFile_MDCparT0_Wir"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_MDCparT0_Wir");
    map_ParamFiles.insert({"mdc_parT0wir_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"mdc_parT0wir_file","./calib/mdc_driftparam/MDC_T0Param_Wir.csv"});

if(Config.IsAvailable("CalibFile_Fiberoffset"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_Fiberoffset");
    map_ParamFiles.insert({"fiber_offset_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"fiber_offset_file","./calib/fiber_offset/fiber_offset.csv"});

if(Config.IsAvailable("CalibFile_FiberTimeoffset"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_FiberTimeoffset");
    map_ParamFiles.insert({"fiber_timeoffset_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"fiber_timeoffset_file","./calib/fiber_offset/fiber_timeoffset.csv"});

if(Config.IsAvailable("CalibFile_FiberAngleoffset"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_FiberAngleoffset");
    map_ParamFiles.insert({"fiber_angleoffset_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"fiber_angleoffset_file","./calib/fiber_offset/fiber_angleoffset.csv"});

if(Config.IsAvailable("Fiber_MFTCor_name"))
  {
    auto tmpStr = Config.Get<std::string>("Fiber_MFTCor_name");
    map_ParamFiles.insert({"fiber_mftcor_file",tmpStr});
  }
  else
    map_ParamFiles.insert({"fiber_mftcor_file","./calib/fiber_mftcor/fiber_mftcor.csv"});

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

if(Config.IsAvailable("CalibFile_MWDCDtDx"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_MWDCDtDx");
    map_ParamFiles.insert({"mwdc_name_dtdx",tmpStr});
  }
  else
    map_ParamFiles.insert({"mwdc_name_dtdx","./calib/mwdc_dtdx/mwdc_dtdx_param.txt"});

if(Config.IsAvailable("CalibFile_MWDCDtDxtable"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_MWDCDtDxtable");
    map_ParamFiles.insert({"mwdc_name_dtdxtable",tmpStr});
  }
  else
    map_ParamFiles.insert({"mwdc_name_dtdxtable","./calib/mwdc_dtdx/mwdc_dtdxtable.root"});
/*
if(Config.IsAvailable("CalibFile_Optics"))
  {
    auto tmpStr = Config.Get<std::string>("CalibFile_Optics");
    map_ParamFiles.insert({"optics_name",tmpStr});
  }
  else
    map_ParamFiles.insert({"optics_name","./calib/optics/optics_par.csv"});
*/

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
  auto vol = gGeoManager->GetVolume("HypHI_RPC_l_log");
  if(vol != nullptr)
    vol->SetMedium(vac);
  vol = gGeoManager->GetVolume("HypHI_RPC_h_log");
  if(vol!=nullptr)
    vol->SetMedium(vac);
  vol = gGeoManager->GetVolume("FMF2_log");
  if(vol!=nullptr)
    vol->SetMedium(vac);
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

  if(InputPar.simParameters != nullptr)
    {
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
    }
  else
    {
      Wasa_FieldMap = true;
      Wasa_FieldMapName = "./field/MagField_default.dat";
      Field_Strength   = Config.IsAvailable("Field_Strength")   ? Config.Get<double>("Field_Strength")   : -1.;
      Target_PositionX = Config.IsAvailable("Target_PositionX") ? Config.Get<double>("Target_PositionX") : 0.;
      Target_PositionY = Config.IsAvailable("Target_PositionY") ? Config.Get<double>("Target_PositionY") : 0.;
      Target_PositionZ = Config.IsAvailable("Target_PositionZ") ? Config.Get<double>("Target_PositionZ") : 196.12;
      Target_Size      = Config.IsAvailable("Target_Size")      ? Config.Get<double>("Target_Size")      : 3.;
    }

  TObjArray* L_vol = gGeoManager->GetListOfVolumes();
  int n_volume     = L_vol->GetEntries();
  for(int n_v = 0; n_v < n_volume; ++n_v)
    {
      std::string name_vol_temp(L_vol->At(n_v)->GetName());
      name_GeoVolumes.push_back(name_vol_temp);
    }
  name_GeoVolumes.push_back("Total");

/*
  std::string upk_setupfiber = InputPar.simexpMetadata->Unpack_SetupFiber;
  Unpack_map_ParamFiles["upk_setupfiber"] = upk_setupfiber;

  std::string upk_channelmapt0 = InputPar.simexpMetadata->Unpack_ChannelMapT0;
  Unpack_map_ParamFiles["upk_channelmapt0"] = upk_channelmapt0;

  std::string upk_channelmappsfe = InputPar.simexpMetadata->Unpack_ChannelMapPSFE;
  Unpack_map_ParamFiles["upk_channelmappsfe"] = upk_channelmappsfe;

  std::string upk_channelmappsbe = InputPar.simexpMetadata->Unpack_ChannelMapPSBE;
  Unpack_map_ParamFiles["upk_channelmappsbe"] = upk_channelmappsbe;

  std::string upk_setuppsb = InputPar.simexpMetadata->Unpack_SetupPSB;
  Unpack_map_ParamFiles["upk_setuppsb"] = upk_setuppsb;

  std::string upk_dtdxtablemwdc = InputPar.simexpMetadata->Unpack_DtDxTableMWDC;
  Unpack_map_ParamFiles["upk_dtdxtablemwdc"] = upk_dtdxtablemwdc;

  std::string upk_celloffsets4wfd1 = InputPar.simexpMetadata->Unpack_CellOffsetS4WFD1;
  Unpack_map_ParamFiles["upk_celloffsets4wfd1"] = upk_celloffsets4wfd1;

  std::string upk_celloffsets2wfd1 = InputPar.simexpMetadata->Unpack_CellOffsetS2WFD1;
  Unpack_map_ParamFiles["upk_celloffsets2wfd1"] = upk_celloffsets2wfd1;

  std::string upk_celloffsets2wfd2 = InputPar.simexpMetadata->Unpack_CellOffsetS2WFD2;
  Unpack_map_ParamFiles["upk_celloffsets2wfd2"] = upk_celloffsets2wfd2;

  std::string upk_celloffsets2wfd3 = InputPar.simexpMetadata->Unpack_CellOffsetS2WFD3;
  Unpack_map_ParamFiles["upk_celloffsets2wfd3"] = upk_celloffsets2wfd3;

  std::string upk_celloffsets2wfd4 = InputPar.simexpMetadata->Unpack_CellOffsetS2WFD4;
  Unpack_map_ParamFiles["upk_celloffsets2wfd4"] = upk_celloffsets2wfd4;

  std::string upk_celloffsets2wfd5 = InputPar.simexpMetadata->Unpack_CellOffsetS2WFD5;
  Unpack_map_ParamFiles["upk_celloffsets2wfd5"] = upk_celloffsets2wfd5;

  std::string upk_channelmapmdc = InputPar.simexpMetadata->Unpack_ChannelMapMDC;
  Unpack_map_ParamFiles["upk_channelmapmdc"] = upk_channelmapmdc;

  std::string upk_physicalmapmdc = InputPar.simexpMetadata->Unpack_PhysicalMapMDC;
  Unpack_map_ParamFiles["upk_physicalmapmdc"] = upk_physicalmapmdc;

  std::string upk_driftparammdc = InputPar.simexpMetadata->Unpack_PhysicalMapMDC;
  Unpack_map_ParamFiles["upk_driftparammdc"] = upk_driftparammdc;
*/

  return 0;
}

int THyphiAttributes::Reload_Para()
{
  std::string Hash(InputPar.simexpMetadata->Hash);

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
  if(Config.IsAvailable("Task_PrimaryVtx_Si"))
    Task_PrimaryVtx_Si = Config.Get<bool>("Task_PrimaryVtx_Si");
  if(Config.IsAvailable("Task_FlatMCOutputML"))
    Task_FlatMCOutputML = Config.Get<bool>("Task_FlatMCOutputML");
  if(Config.IsAvailable("Task_CheckFiberXUV"))
    Task_CheckFiberXUV = Config.Get<bool>("Task_CheckFiberXUV");
  if(Config.IsAvailable("Task_CheckFiberTrack"))
    Task_CheckFiberTrack = Config.Get<bool>("Task_CheckFiberTrack");
  if(Config.IsAvailable("Task_FragmentFinder"))
    Task_FragmentFinder = Config.Get<bool>("Task_FragmentFinder");
  if(Config.IsAvailable("Task_WASAFinder"))
    Task_WASAFinder = Config.Get<bool>("Task_WASAFinder");
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
  if(Config.IsAvailable("Task_DecayVtx_pi+"))
    Task_DecayVtx_piplus = Config.Get<bool>("Task_DecayVtx_pi+");
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
    if(s == "Task_PrimaryVtx_Si")
      Task_Order.push_back(TASKPRIMARYVTX_SI);
	  if(s == "Task_FlatMCOutputML")
	    Task_Order.push_back(TASKFLATMCOUTPUTML);
    if(s == "Task_CheckFiberXUV")
      Task_Order.push_back(TASKCHECKFIBERXUV);
	  if(s == "Task_CheckFiberTrack")
	    Task_Order.push_back(TASKCHECKFIBERTRACK);
    if(s == "Task_FragmentFinder")
      Task_Order.push_back(TASKFRAGMENTFINDER);
    if(s == "Task_WASAFinder")
      Task_Order.push_back(TASKWASAFINDER);
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
    if(s == "Task_DecayVtx_pi+")
      Task_Order.push_back(TASKDECAYVTX_PIPLUS);
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
  out.RunTask.Task_PrimaryVtx_Si   = TaskConfig.Task_PrimaryVtx_Si;
  out.RunTask.Task_FlatMCOutputML  = TaskConfig.Task_FlatMCOutputML;
  out.RunTask.Task_FragmentFinder  = TaskConfig.Task_FragmentFinder;
  out.RunTask.Task_WASAFinder      = TaskConfig.Task_WASAFinder;
  out.RunTask.Task_BayesFinder     = TaskConfig.Task_BayesFinder;
  out.RunTask.Task_RiemannFinder   = TaskConfig.Task_RiemannFinder;
  out.RunTask.Task_FinderCM        = TaskConfig.Task_FinderCM;
  out.RunTask.Task_FindingPerf     = TaskConfig.Task_FindingPerf;
  out.RunTask.Task_CheckRZ         = TaskConfig.Task_CheckRZ;
  out.RunTask.Task_KalmanDAF       = TaskConfig.Task_KalmanDAF;
  out.RunTask.Task_DecayVtx        = TaskConfig.Task_DecayVtx;
  out.RunTask.Task_DecayVtx_piplus = TaskConfig.Task_DecayVtx_piplus;
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
  out.RunTaskAttr.PV_RealXUVComb     = PV_RealXUVComb;
  out.RunTaskAttr.PV_RealPrimTrack   = PV_RealPrimTrack;
  out.RunTaskAttr.RZ_ChangeMiniFiber = RZ_ChangeMiniFiber;
  out.RunTaskAttr.RZ_MDCProlate      = RZ_MDCProlate;
  out.RunTaskAttr.RZ_MDCWire2        = RZ_MDCWire2;
  out.RunTaskAttr.RZ_MDCBiasCorr     = RZ_MDCBiasCorr;
  out.RunTaskAttr.WASAFinder_perfect = WASAFinder_perfect;
  out.RunTaskAttr.WASAFinder_PSBHits = WASAFinder_PSBHits;
  out.RunTaskAttr.WASAFinder_PSFEHits= WASAFinder_PSFEHits;
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
