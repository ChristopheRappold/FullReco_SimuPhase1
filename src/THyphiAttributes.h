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

#include "spdlog/logger.h"

struct DataSim
{
  std::vector<std::string>* nameDet;
  std::map<std::string, double>* simParameters;
};

class THyphiAttributes
{

  public:
  int Nb_CPU;
  int Nb_Fraction;
  int NEvent;

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

  std::vector<std::string> name_GeoVolumes;

  bool beam_only;
  bool Debug_DAF;
  bool DoNoMaterial;

  bool Task_CheckField;
  bool Task_PrimaryVtx;
  bool Task_FlatMCOutputML;
  bool Task_BayesFinder;
  bool Task_RiemannFinder;
  bool Task_FinderCM;
  bool Task_CheckRZ;
  bool Task_KalmanDAF;
  bool Task_DecayVtx;

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

  int Init_Para();
};

#endif
