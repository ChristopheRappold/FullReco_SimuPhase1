#ifndef FullRecoTaskMT_h
#define FullRecoTaskMT_h

#include <list>
#include <map>
#include <string>

#include "Ana_Event/Ana_EventNew_v16.hh"
#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Hist.hh"
#include "Debug.hh"
#include "EventG4Sol/TG4Sol_Event.hh"
#include "EventG4Sol/TG4Sol_Hit.hh"
#include "FullRecoEvent.hh"
#include "TDataBuilder.h"
#include "TDataProcess.h"
#include "TDataMerger.h"
#include "THyphiAttributes.h"

#include "TTree.h"

#ifdef ROOT6
#include "TTreeReaderArray.h"
#endif

class FullRecoEvent;
class Ana_Hist;

class FullRecoTaskMT
{

public:
  FullRecoTaskMT() = delete;
  FullRecoTaskMT(const FullRecoConfig& conf, const DataSim& In);
  ~FullRecoTaskMT();
  
  int Run(Long64_t startEvent, Long64_t stopEvent, TTree* InTree, const TG4Sol_Event& ev,
	  const std::vector<TClonesArray*>& hits, TTree* OutTree, MCAnaEventG4Sol* OutEvent);
  void AttachHisto(Ana_Hist* h);

private:
  const THyphiAttributes Attributes;

  std::vector<FullRecoEvent> REvent;
  Ana_Hist* AnaHisto;

  std::vector<TDataBuilder*> list_det_build;
  std::vector<TDataMerger*> list_merger_out;

  std::vector<TDataProcess<FullRecoEvent, int>*> list_processMC_MT;

  int EventProcess(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutEvent);
};

#endif
