#ifndef THYPHIATTRIBUTE
#define THYPHIATTRIBUTE
//#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>
//#include "TRandom3.h"
#include "FullRecoConfig.hh"

#include "FairBase/FrsSksHypFieldMapFull.h"
#include "FairBase/HypConstField.h"
#include "FairBase/HypFieldMap.h"
#include "FairBase/HypFieldMapFull.h"

#include "FairBase/FrsSolenoidHypField.h"
#include "FairBase/WasaSolenoidFieldMap.h"

#include "Ana_Event/AnaEvent_Metadata.hh"

#include "spdlog/logger.h"

#include "sqlite_orm/sqlite_orm.h"

struct DataSim
{
  std::vector<std::string>* nameDet;
  std::map<std::string, double>* simParameters;

  AnaEvent_Metadata* previousMeta;
};

struct MT
{
  enum MTtype : int
    {
      SingleTh = 0,
      MultiTh,
      ZeroMQ,
    };
  MTtype RunType;

  bool IsMain;
  int NQueue;
  int NBuilder;
  int NKalman;
  int NMerger;

  std::string addr_initEvent;

  std::string addr_frontBuilder;
  std::string addr_backFitter;

  std::string addr_frontFitter;
  std::string addr_backMerger;

  std::string addr_frontMerger;
  std::string addr_backEnd;

  std::string addr_control;
  std::string addr_monitor;

  void Init(const FullRecoConfig& Config);

};

struct Task
{
  bool Task_ReStart = false;
  bool Task_CheckField = true;
  bool Task_PrimaryVtx = false;
  bool Task_FlatMCOutputML = false;
  bool Task_BayesFinder = false;
  bool Task_RiemannFinder = false;
  bool Task_FinderCM = false;
  bool Task_FindingPerf = false;
  bool Task_CheckRZ = true;
  bool Task_KalmanDAF = true;
  bool Task_DecayVtx = false;

  void Init(const FullRecoConfig& Config);
};

// struct KF
// {


// };

struct AttrOut
{
  std::string Hash;

  std::string NameIn;
  std::string NameOut;
  std::string DateOfRun;

  int Nb_CPU;
  int Nb_Fraction;
  int NEvent;

  int MT_RunType;
  bool MT_IsMain;
  int MT_NQueue;
  int MT_NBuilder;
  int MT_NKalman;
  int MT_NMerger;
  std::string MT_addr_initEvent;
  std::string MT_addr_frontBuilder;
  std::string MT_addr_backFitter;
  std::string MT_addr_frontFitter;
  std::string MT_addr_backMerger;
  std::string MT_addr_frontMerger;
  std::string MT_addr_backEnd;
  std::string MT_addr_control;
  std::string MT_addr_monitor;

  bool Task_CheckField;
  bool Task_PrimaryVtx;
  bool Task_FlatMCOutputML;
  bool Task_BayesFinder;
  bool Task_RiemannFinder;
  bool Task_FinderCM;
  bool Task_FindingPerf;
  bool Task_CheckRZ;
  bool Task_KalmanDAF;
  bool Task_DecayVtx;
  bool Task_ReStart;

  bool G4_simu;
  bool G4_TimeResolution;
  bool G4_GeoResolution;

  bool back_tracking;

  double Target_PositionX;
  double Target_PositionY;
  double Target_PositionZ;

  double Target_Size;
  double Field_Strength;

  int Wasa_Side;
  bool Wasa_FieldMap;
  std::string Wasa_FieldMapName;

  std::vector<char> name_GeoVolumes;

  bool beam_only;
  bool Debug_DAF;
  bool DoNoMaterial;

  bool RZ_ChangeMiniFiber;
  bool RZ_MDCProlate;
  bool RZ_MDCWire2;
  bool RZ_MDCBiasCorr;

  bool KF_Kalman;
  bool KF_KalmanSqrt;
  bool KF_KalmanRef;
  bool KF_DAFRef;
  bool KF_DAF;

  int KF_NbCentralCut;
  int KF_NbMiniFiberCut;

  std::string FlatML_namefile;
  std::string DataML_Out;

  bool RF_OutputEvents;

  std::vector<char> ConfigJson;
};

inline auto InitStorage()
{
  using namespace sqlite_orm;
  return make_storage(
		      "WasaFRSexp.sqlite",
		      make_table(
				 "AttrRunDone", make_column("HashId", &AttrOut::Hash), make_column("NameIn", &AttrOut::NameIn),
				 make_column("NameOut", &AttrOut::NameOut), make_column("DateOfRun", &AttrOut::DateOfRun),
				 make_column("Nb_CPU", &AttrOut::Nb_CPU), make_column("Nb_Fraction", &AttrOut::Nb_Fraction),
				 make_column("NEvent", &AttrOut::NEvent), make_column("MT_RunType", &AttrOut::MT_RunType),
				 make_column("MT_IsMain", &AttrOut::MT_IsMain), make_column("MT_NQueue", &AttrOut::MT_NQueue),
				 make_column("MT_NBuilder", &AttrOut::MT_NBuilder), make_column("MT_NKalman", &AttrOut::MT_NKalman),
				 make_column("MT_NMerger", &AttrOut::MT_NMerger),
				 make_column("MT_add_initEvent", &AttrOut::MT_addr_initEvent),
				 make_column("MT_add_frontBuilder", &AttrOut::MT_addr_frontBuilder),
				 make_column("MT_add_backFitter", &AttrOut::MT_addr_backFitter),
				 make_column("MT_add_frontFitter", &AttrOut::MT_addr_frontFitter),
				 make_column("MT_add_backMerger", &AttrOut::MT_addr_backMerger),
				 make_column("MT_add_frontMerger", &AttrOut::MT_addr_frontMerger),
				 make_column("MT_add_backEnd", &AttrOut::MT_addr_backEnd),
				 make_column("MT_add_control", &AttrOut::MT_addr_control),
				 make_column("MT_add_monitor", &AttrOut::MT_addr_monitor),
				 make_column("Task_CheckField", &AttrOut::Task_CheckField),
				 make_column("Task_PrimaryVtx", &AttrOut::Task_PrimaryVtx),
				 make_column("Task_FlatMCOutputML", &AttrOut::Task_FlatMCOutputML),
				 make_column("Task_BayesFinder", &AttrOut::Task_BayesFinder),
				 make_column("Task_RiemannFinder", &AttrOut::Task_RiemannFinder),
				 make_column("Task_FinderCM", &AttrOut::Task_FinderCM),
				 make_column("Task_FindingPerf", &AttrOut::Task_FindingPerf),
				 make_column("Task_CheckRZ", &AttrOut::Task_CheckRZ),
				 make_column("Task_KalmanDAF", &AttrOut::Task_KalmanDAF),
				 make_column("Task_DecayVtx", &AttrOut::Task_DecayVtx),
				 make_column("Task_ReStart", &AttrOut::Task_ReStart), make_column("G4_simu", &AttrOut::G4_simu),
				 make_column("G4_TimeResolution", &AttrOut::G4_TimeResolution),
				 make_column("G4_GeoResolution", &AttrOut::G4_GeoResolution),
				 make_column("back_tracking", &AttrOut::back_tracking),
				 make_column("Target_PositionX", &AttrOut::Target_PositionX),
				 make_column("Target_PositionY", &AttrOut::Target_PositionY),
				 make_column("Target_PositionZ", &AttrOut::Target_PositionZ),
				 make_column("Target_Size", &AttrOut::Target_Size),
				 make_column("Field_Strength", &AttrOut::Field_Strength),
				 make_column("Wasa_Side", &AttrOut::Wasa_Side), make_column("Wasa_FieldMap", &AttrOut::Wasa_FieldMap),
				 make_column("Wasa_FieldMapName", &AttrOut::Wasa_FieldMapName),
				 make_column("name_GeoVolumes", &AttrOut::name_GeoVolumes),
				 make_column("beam_only", &AttrOut::beam_only), make_column("Debug_DAF", &AttrOut::Debug_DAF),
				 make_column("DoNoMaterial", &AttrOut::DoNoMaterial),
				 make_column("RZ_ChangeMiniFiber", &AttrOut::RZ_ChangeMiniFiber),
				 make_column("RZ_MDCProlate", &AttrOut::RZ_MDCProlate),
				 make_column("RZ_MDCWire2", &AttrOut::RZ_MDCWire2),
				 make_column("RZ_MDCBiasCorr", &AttrOut::RZ_MDCBiasCorr),
				 make_column("KF_Kalman", &AttrOut::KF_Kalman), make_column("KF_KalmanSqrt", &AttrOut::KF_KalmanSqrt),
				 make_column("KF_KalmanRef", &AttrOut::KF_KalmanRef), make_column("KF_DAFRef", &AttrOut::KF_DAFRef),
				 make_column("KF_DAF", &AttrOut::KF_DAF), make_column("KF_NbCentralCut", &AttrOut::KF_NbCentralCut),
				 make_column("KF_MbMiniFiberCut", &AttrOut::KF_NbMiniFiberCut),
				 make_column("FlatML_namefile", &AttrOut::FlatML_namefile),
				 make_column("DataML_Out", &AttrOut::DataML_Out),
				 make_column("RF_OutputEvents", &AttrOut::RF_OutputEvents),
				 make_column("ConfigJson", &AttrOut::ConfigJson)));

}



class THyphiAttributes
{

  public:

  std::string NameIn;
  std::string NameOut;
  std::string DateOfRun;

  std::string Hash;

  int Nb_CPU;
  int Nb_Fraction;
  int NEvent;

  MT MTsetting;

  Task TaskConfig;

  bool G4_simu;
  bool G4_TimeResolution;
  bool G4_GeoResolution;

  bool back_tracking;

  double Target_PositionX;
  double Target_PositionY;
  double Target_PositionZ;

  double Target_Size;
  double Field_Strength;

  int Wasa_Side;
  bool Wasa_FieldMap;
  std::string Wasa_FieldMapName;

  std::vector<std::string> name_GeoVolumes;

  bool beam_only;
  bool Debug_DAF;
  bool DoNoMaterial;

  bool RZ_ChangeMiniFiber;
  bool RZ_MDCProlate;
  bool RZ_MDCWire2;
  bool RZ_MDCBiasCorr;

  bool KF_Kalman;
  bool KF_KalmanSqrt;
  bool KF_KalmanRef;
  bool KF_DAFRef;
  bool KF_DAF;

  int KF_NbCentralCut;
  int KF_NbMiniFiberCut;

  std::string FlatML_namefile;
  std::string DataML_Out;

  bool RF_OutputEvents;

  FairField* Field;

  const FullRecoConfig& Config;
  const DataSim& InputPar;

  std::shared_ptr<spdlog::logger> _logger;
  
  THyphiAttributes() = delete;
  THyphiAttributes(const FullRecoConfig& conf, const DataSim& InputParameters);
  ~THyphiAttributes() = default;

  void SetOut(AttrOut& out) const;

  int Init_Para();
  int Reload_Para();
  int CrossCheckConfig();

  void SaveToDatabase() const;
};

#endif
