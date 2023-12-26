#ifndef THYPHIATTRIBUTE
#define THYPHIATTRIBUTE
//#include <iostream>
#include <list>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include <fstream>
//#include "TRandom3.h"

#include "TFile.h"
#include "TCutG.h"

#include "FullRecoConfig.hh"

#include "FairBase/FrsSksHypFieldMapFull.h"
#include "FairBase/HypConstField.h"
#include "FairBase/HypFieldMap.h"
#include "FairBase/HypFieldMapFull.h"

#include "FairBase/FrsSolenoidHypField.h"
#include "FairBase/WasaSolenoidFieldMap.h"

#include "Ana_Event/AnaEvent_Metadata.hh"
//#include "EventWASAUnpack/TFRSParameter.h"

#include "HitAna/ParaManager.hh"

#include "spdlog/logger.h"

#include "sqlite_orm/sqlite_orm.h"

struct DataSimExp
{
  std::vector<std::string>* nameDet;  
  std::map<std::string, double>* simParameters;

  AnaEvent_Metadata* simexpMetadata;
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
  bool Task_CheckField = false;
  bool Task_FiberHitFinder = true;
  bool Task_PrimaryVtx = true;
  bool Task_PrimaryVtx_Si = false;
  bool Task_FlatMCOutputML = false;
  bool Task_CheckFiberXUV = false;
  bool Task_CheckFiberTrack = false;
  bool Task_FragmentFinder = false;
  bool Task_WASAFinder = false;
  bool Task_GNNFinder = false;
  bool Task_BayesFinder = false;
  bool Task_RiemannFinder = false;
  bool Task_FinderCM = false;
  bool Task_FindingPerf = false;
  bool Task_CheckRZ = true;
  bool Task_RPhiZTrackMDC = false;
  bool Task_KalmanDAF = true;
  bool Task_KalmanDAFPID = true;
  bool Task_DecayVtx = false;
  bool Task_DecayVtx_piplus = false;
  enum Task_Id
  {
    TASKRESTART = 0, TASKCHECKFIELD, TASKFIBERHITFINDER, TASKPRIMARYVTX, TASKPRIMARYVTX_SI, TASKFLATMCOUTPUTML, TASKCHECKFIBERXUV, TASKCHECKFIBERTRACK, TASKFRAGMENTFINDER, TASKWASAFINDER, TASKGNNFINDER, TASKBAYESFINDER, TASKRIEMANNFINDER, TASKFINDERCM, TASKFINDINGPERF, TASKCHECKRZ, TASKRPHIZTRACKMDC,
    TASKKALMANDAF, TASKKALMANDAFPID, TASKDECAYVTX, TASKDECAYVTX_PIPLUS, NBTASKID
  };

  std::vector<Task_Id> Task_Order = {TASKCHECKFIELD, TASKFIBERHITFINDER, TASKPRIMARYVTX, TASKPRIMARYVTX_SI, TASKCHECKFIBERXUV, TASKCHECKFIBERTRACK, TASKFRAGMENTFINDER, TASKWASAFINDER, TASKGNNFINDER, TASKBAYESFINDER, TASKRIEMANNFINDER, TASKFINDINGPERF, TASKCHECKRZ, TASKRPHIZTRACKMDC, TASKKALMANDAF, TASKKALMANDAFPID, TASKDECAYVTX, TASKDECAYVTX_PIPLUS, TASKFLATMCOUTPUTML};

  void Init(const FullRecoConfig& Config);
};


struct RunDoneDef
  {
    std::string Hash;

    std::string NameIn;
    std::string NameOut;
    std::string DateOfRun;

    int Nb_CPU;
    int Nb_Fraction;
    int NEvent;
  };

struct RunMultiDef
  {
    std::string Hash;

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
};
struct RunTaskDef
{
  std::string Hash;

  bool Task_CheckField;
  bool Task_FiberHitFinder;
  bool Task_PrimaryVtx;
  bool Task_PrimaryVtx_Si;
  bool Task_FlatMCOutputML;
  bool Task_CheckFiberXUV;
  bool Task_CheckFiberTrack;
  bool Task_FragmentFinder;
  bool Task_WASAFinder;
  bool Task_GNNFinder;
  bool Task_BayesFinder;
  bool Task_RiemannFinder;
  bool Task_FinderCM;
  bool Task_FindingPerf;
  bool Task_CheckRZ;
  bool Task_RPhiZTrackMDC;
  bool Task_KalmanDAF;
  bool Task_KalmanDAFPID;
  bool Task_DecayVtx;
  bool Task_DecayVtx_piplus;
  bool Task_ReStart;
};

struct RunAttrGeneralDef
{
  std::string Hash;

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
};

struct RunTaskAttrDef
{
  std::string Hash;

  bool Debug_DAF;
  bool DoNoMaterial;

  bool PV_RealXUVComb;
  bool PV_RealPrimTrack;
  bool CFT_RZfit;

  bool RZ_ChangeMiniFiber;
  bool RZ_MDCProlate;
  bool RZ_MDCWire2;
  bool RZ_MDCBiasCorr;

  bool WF_perfect;
  bool WF_PSBHits;
  bool WF_PSFEHits;

  bool GNN_Text;
  std::string GNN_Node;

  bool KF_Kalman;
  bool KF_KalmanSqrt;
  bool KF_KalmanRef;
  bool KF_DAFRef;
  bool KF_DAF;
  bool KF_G4e;

  std::string G4e_FullProp;
  std::string G4e_Basf2List;
  std::string G4e_ExactJac;
  std::string G4e_MaxEnergyLoss;
  
  int KF_NbCentralCut;
  int KF_NbMiniFiberCut;

  std::string FlatML_namefile;
  std::string DataML_Out;

  bool RF_OutputEvents;
  bool CFT_OutputEvents;
};

struct RunFullConfigDef
{
  std::string Hash;
  std::vector<char> ConfigJson;
};

struct AttrOut
{
  RunDoneDef RunDone;
  RunMultiDef RunMulti;
  RunTaskDef RunTask;
  RunAttrGeneralDef RunAttrGeneral;
  RunTaskAttrDef RunTaskAttr;
  RunFullConfigDef RunFullConfig;
};

inline auto InitStorage()
{
  using namespace sqlite_orm;
  return make_storage(
		      "WasaFRSexp.sqlite",
		      make_table(
				 "RunDone",
				 make_column("HashId", &RunDoneDef::Hash),
				 make_column("NameIn", &RunDoneDef::NameIn),
				 make_column("NameOut", &RunDoneDef::NameOut),
				 make_column("DateOfRun", &RunDoneDef::DateOfRun),
				 make_column("Nb_CPU", &RunDoneDef::Nb_CPU),
				 make_column("Nb_Fraction", &RunDoneDef::Nb_Fraction),
				 make_column("NEvent", &RunDoneDef::NEvent)
				 ),
		      make_table(
				 "RunMultiThread",
				 make_column("HashId", &RunMultiDef::Hash),
				 make_column("MT_RunType", &RunMultiDef::MT_RunType),
				 make_column("MT_IsMain", &RunMultiDef::MT_IsMain),
				 make_column("MT_NQueue", &RunMultiDef::MT_NQueue),
				 make_column("MT_NBuilder", &RunMultiDef::MT_NBuilder),
				 make_column("MT_NKalman", &RunMultiDef::MT_NKalman),
				 make_column("MT_NMerger", &RunMultiDef::MT_NMerger),
				 make_column("MT_add_initEvent", &RunMultiDef::MT_addr_initEvent),
				 make_column("MT_add_frontBuilder", &RunMultiDef::MT_addr_frontBuilder),
				 make_column("MT_add_backFitter", &RunMultiDef::MT_addr_backFitter),
				 make_column("MT_add_frontFitter", &RunMultiDef::MT_addr_frontFitter),
				 make_column("MT_add_backMerger", &RunMultiDef::MT_addr_backMerger),
				 make_column("MT_add_frontMerger", &RunMultiDef::MT_addr_frontMerger),
				 make_column("MT_add_backEnd", &RunMultiDef::MT_addr_backEnd),
				 make_column("MT_add_control", &RunMultiDef::MT_addr_control),
				 make_column("MT_add_monitor", &RunMultiDef::MT_addr_monitor)
				 ),
		      make_table(
				 "RunTask",
				 make_column("HashId", &RunTaskDef::Hash),
				 make_column("Task_CheckField", &RunTaskDef::Task_CheckField),
				 make_column("Task_PrimaryVtx", &RunTaskDef::Task_PrimaryVtx),
				 make_column("Task_PrimaryVtx_Si", &RunTaskDef::Task_PrimaryVtx_Si),
				 make_column("Task_FlatMCOutputML", &RunTaskDef::Task_FlatMCOutputML),
				 make_column("Task_FragmentFinder", &RunTaskDef::Task_FragmentFinder),
				 make_column("Task_WASAFinder", &RunTaskDef::Task_WASAFinder),
         make_column("Task_GNNFinder", &RunTaskDef::Task_GNNFinder),
				 make_column("Task_BayesFinder", &RunTaskDef::Task_BayesFinder),
				 make_column("Task_RiemannFinder", &RunTaskDef::Task_RiemannFinder),
				 make_column("Task_FinderCM", &RunTaskDef::Task_FinderCM),
				 make_column("Task_FindingPerf", &RunTaskDef::Task_FindingPerf),
				 make_column("Task_CheckRZ", &RunTaskDef::Task_CheckRZ),
				 make_column("Task_RPhiZTrackMDC", &RunTaskDef::Task_RPhiZTrackMDC),
				 make_column("Task_KalmanDAF", &RunTaskDef::Task_KalmanDAF),
				 make_column("Task_KalmanDAFPID", &RunTaskDef::Task_KalmanDAFPID),
				 make_column("Task_DecayVtx", &RunTaskDef::Task_DecayVtx),
         make_column("Task_DecayVtx_piPlus", &RunTaskDef::Task_DecayVtx_piplus),
				 make_column("Task_ReStart", &RunTaskDef::Task_ReStart)
				 ),
		      make_table(
				 "RunAttrGeneral",
				 make_column("HashId", &RunAttrGeneralDef::Hash),
				 make_column("G4_simu", &RunAttrGeneralDef::G4_simu),
				 make_column("G4_TimeResolution", &RunAttrGeneralDef::G4_TimeResolution),
				 make_column("G4_GeoResolution", &RunAttrGeneralDef::G4_GeoResolution),
				 make_column("back_tracking", &RunAttrGeneralDef::back_tracking),
				 make_column("Target_PositionX", &RunAttrGeneralDef::Target_PositionX),
				 make_column("Target_PositionY", &RunAttrGeneralDef::Target_PositionY),
				 make_column("Target_PositionZ", &RunAttrGeneralDef::Target_PositionZ),
				 make_column("Target_Size", &RunAttrGeneralDef::Target_Size),
				 make_column("Field_Strength", &RunAttrGeneralDef::Field_Strength),
				 make_column("Wasa_Side", &RunAttrGeneralDef::Wasa_Side),
				 make_column("Wasa_FieldMap", &RunAttrGeneralDef::Wasa_FieldMap),
				 make_column("Wasa_FieldMapName", &RunAttrGeneralDef::Wasa_FieldMapName),
				 make_column("name_GeoVolumes", &RunAttrGeneralDef::name_GeoVolumes),
				 make_column("beam_only", &RunAttrGeneralDef::beam_only)
				 ),
		      make_table(
				 "RunTaskAttr",
				 make_column("HashId", &RunTaskAttrDef::Hash),
				 make_column("Debug_DAF", &RunTaskAttrDef::Debug_DAF),
				 make_column("DoNoMaterial", &RunTaskAttrDef::DoNoMaterial),
				 make_column("PV_RealXUVComb", &RunTaskAttrDef::PV_RealXUVComb),
				 make_column("PV_RealPrimTrack", &RunTaskAttrDef::PV_RealPrimTrack),
				 make_column("CFT_RZfit", &RunTaskAttrDef::CFT_RZfit),
				 make_column("RZ_ChangeMiniFiber", &RunTaskAttrDef::RZ_ChangeMiniFiber),
				 make_column("RZ_MDCProlate", &RunTaskAttrDef::RZ_MDCProlate),
				 make_column("RZ_MDCWire2", &RunTaskAttrDef::RZ_MDCWire2),
				 make_column("RZ_MDCBiasCorr", &RunTaskAttrDef::RZ_MDCBiasCorr),
				 make_column("WF_perfect", &RunTaskAttrDef::WF_perfect),
				 make_column("WF_PSBHits", &RunTaskAttrDef::WF_PSBHits),
				 make_column("WF_PSFEHits", &RunTaskAttrDef::WF_PSFEHits),
				 make_column("GNN_Text", &RunTaskAttrDef::GNN_Text),
				 make_column("GNN_Node", &RunTaskAttrDef::GNN_Node),
				 make_column("KF_Kalman", &RunTaskAttrDef::KF_Kalman),
				 make_column("KF_KalmanSqrt", &RunTaskAttrDef::KF_KalmanSqrt),
				 make_column("KF_KalmanRef", &RunTaskAttrDef::KF_KalmanRef),
				 make_column("KF_DAFRef", &RunTaskAttrDef::KF_DAFRef),
				 make_column("KF_DAF", &RunTaskAttrDef::KF_DAF),
				 make_column("KF_NbCentralCut", &RunTaskAttrDef::KF_NbCentralCut),
				 make_column("KF_MbMiniFiberCut", &RunTaskAttrDef::KF_NbMiniFiberCut),
				 make_column("FlatML_namefile", &RunTaskAttrDef::FlatML_namefile),
				 make_column("DataML_Out", &RunTaskAttrDef::DataML_Out),
				 make_column("RF_OutputEvents", &RunTaskAttrDef::RF_OutputEvents),
				 make_column("CFT_OutputEvents", &RunTaskAttrDef::CFT_OutputEvents)
				 ),
		      make_table(
				 "RunFullConfig",
				 make_column("HashId", &RunFullConfigDef::Hash),
				 make_column("ConfigJson", &RunFullConfigDef::ConfigJson)
				 )
		      );

}



class THyphiAttributes
{

  public:

  std::string NameIn;
  std::string NameOut;
  std::string DateOfRun;

  std::string Hash;

  std::unordered_map<std::string,std::string> map_ParamFiles;
  //std::unordered_map<std::string,std::string> Unpack_map_ParamFiles;

  // Optics parameters
  std::string optics_name;
  std::map<std::string, std::vector<float>> optics_par;
  double optics_s2z;

  // PID Cuts
  std::string pi_cut_name;
  TCutG *cut_pi;
  
  // MWDC Par load
  int mwdc_par_fitorhist; // 0 -> fit ; 1 -> Hist

  // MFT
  double fiber_mft1_pos_z;
  double fiber_mft2_pos_z;

  // PSB
  double psb_pos_x;
  double psb_pos_y;
  double psb_pos_z;
  double psb_rot_z;
  double cut_psb_phi;
  double cut_psb_z;
  double cut_phi_fm;

  bool flag_dup_trackhit;
  bool flag_dup_trackhit_mdc;
  bool flag_trackhit_inclusive;
  
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

  bool PV_RealXUVComb;
  bool PV_RealPrimTrack;
  bool CFT_RZfit;

  bool RZ_ChangeMiniFiber;
  bool RZ_MDCProlate;
  bool RZ_MDCWire2;
  bool RZ_MDCBiasCorr;

  bool WF_perfect;
  bool WF_PSBHits;
  bool WF_PSFEHits;

  bool GNN_Text;
  std::string GNN_Node;

  bool KF_Kalman;
  bool KF_KalmanSqrt;
  bool KF_KalmanRef;
  bool KF_DAFRef;
  bool KF_DAF;
  bool KF_G4e;

  std::string G4e_FullProp;
  std::string G4e_Basf2List;
  std::string G4e_ExactJac;
  std::string G4e_MaxEnergyLoss;
  std::string G4e_Verbose;

  int KF_NbCentralCut;
  int KF_NbMiniFiberCut;

  std::string StudyCase;

  std::string FlatML_namefile;
  std::string DataML_Out;

  bool RF_OutputEvents;
  bool CFT_OutputEvents;

  FairField* Field;

  const FullRecoConfig& Config;
  const DataSimExp& InputPar;

  std::unique_ptr<ParaManager> par;

  std::shared_ptr<spdlog::logger> _logger;
  
  THyphiAttributes() = delete;
  THyphiAttributes(const FullRecoConfig& conf, const DataSimExp& InputParameters);
  ~THyphiAttributes() = default;

  void SetOut(AttrOut& out) const;

  int Init_Para();
  int Reload_Para();
  int CrossCheckConfig();

  void SaveToDatabase() const;
};

#endif
