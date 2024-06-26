#ifndef FullRecoTask_h
#define FullRecoTask_h


#include <map>
#include <string>
#include <list>
//#include <iostream>

#include "EventG4Sol/TG4Sol_Event.hh"
#include "EventG4Sol/TG4Sol_Hit.hh"
#include "Ana_Event/Ana_EventNew_v16.hh"
#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "Ana_Event/AnaEvent_Metadata.hh"
#include "EventWASAUnpack/WASAUnpackBranch.hh"

#include "FullRecoEvent.hh"
#include "TDataProcess.h"
#include "TDataBuilder.h"

#include "FullRecoConfig.hh"
#include "THyphiAttributes.h"
#include "Ana_Hist.hh"
#include "Debug.hh"

#ifdef ROOT6
#include "TTreeReaderArray.h"
#endif

class FullRecoEvent;
class Ana_Hist;
class FullRecoConfig;

namespace recotask {

  template <typename T, typename = int>
  struct HasMC_Particle : std::false_type { };

  template <typename T>
  struct HasMC_Particle <T, decltype((void) T::fMC_Particle, 0)> : std::true_type { };


};


template<class TEOut>
class FullRecoTask
{

public :
  FullRecoTask() = delete;
  FullRecoTask(const FullRecoConfig& config, const DataSimExp& In);

 ~FullRecoTask();


#ifdef ROOT6
  int EventLoop(const TG4Sol_Event& ev, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, TEOut* OutTree);
#else
  int EventLoop(const TG4Sol_Event& ev, const std::vector<TClonesArray*>& hits, TEOut* OutTree);
#endif  

#ifdef ROOT6
  int EventLoop(const TEOut& RestartEvent, TEOut* OutTree);
#else
  int EventLoop(TEOut* RestartEvent, TEOut* OutTree);
#endif

#ifdef ROOT6
  int EventLoop(const EventWASAUnpack& UnpackEvent, TEOut* OutTree);
#else
  int EventLoop(const EventWASAUnpack& UnpackEvent, TEOut* OutTree);
#endif


  void AttachHisto(Ana_Hist* h);
  void SetEventMetadata(AnaEvent_Metadata& metadata);

private :

  const THyphiAttributes Attributes;

  FullRecoEvent REvent;
  Ana_Hist* AnaHisto;
  TDataBuilder* det_build;

  std::list<TDataProcess<FullRecoEvent,TEOut>*> list_processMC;

  int EventProcess(FullRecoEvent& RecoEvent,TEOut* OutTree);
  

};

#endif
