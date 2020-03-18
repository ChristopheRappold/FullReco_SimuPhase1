#ifndef FullRecoTaskZMQ_h
#define FullRecoTaskZMQ_h

#include <list>
#include <map>
#include <string>

#include "Ana_Event/Ana_EventNew_v16.hh"
#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Hist.hh"
#include "Debug.hh"
#include "EventG4Sol/TG4Sol_Event.hh"
#include "EventG4Sol/TG4Sol_Hit.hh"
#include "FullRecoEventZMQ.hh"
#include "TDataBuilder.h"
#include "TDataProcess.h"
#include "TDataMerger.h"
#include "THyphiAttributes.h"

#include "TTree.h"

#ifdef ROOT6
#include "TTreeReaderArray.h"
#endif

#include "ThreadModelZMQ.hh"
#include "ThreadTasks.hh"
#include "TBuildDetectorLayerPlaneDAF_ZMQ.h"
#include "TKalmanFilter_DAF_ZMQ.h"
#include "TMergerOutput_ZMQ.h"


class Ana_Hist;

using ProcessFitter = TDataProcess<ZMQ::DataBuilderOut, ZMQ::DataFitterOut>;

class FullRecoTaskZMQ
{

public:
  FullRecoTaskZMQ() = delete;
  FullRecoTaskZMQ(const FullRecoConfig& conf, const DataSim& In);
  ~FullRecoTaskZMQ();
  
  int Run(Long64_t startEvent, Long64_t stopEvent, TTree* InTree, const TG4Sol_Event& ev,
	  const std::vector<TClonesArray*>& hits, TTree* OutTree, MCAnaEventG4Sol* OutEvent);
  void AttachHisto(Ana_Hist* h);

private:
  const THyphiAttributes Attributes;

  //std::vector<FullRecoEvent> REvent;
  Ana_Hist* AnaHisto;

  std::vector<TDataBuilder*> list_det_build;
  std::vector<TDataMerger*> list_merger_out;

  std::vector<ProcessFitter*> list_processMC_MT;

  bool IsMain;
  int sizeBuilder;
  int sizeKalman;
  int sizeMerger;
  
  std::string addr_initEvent = "";

  std::string addr_frontBuilder = ""; 
  std::string addr_backFitter = "";

  std::string addr_frontFitter = "";
  std::string addr_backMerger = "";

  std::string addr_frontMerger = "";
  std::string addr_backEnd = "";

  std::string addr_control = "";
  std::string addr_monitor = "";
  
  zmq::context_t contextTask;  

  //int EventProcess(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutEvent);
};


#endif
