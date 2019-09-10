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

class FullRecoTask
{

public :
  FullRecoTask() = delete;
  FullRecoTask(const FullRecoConfig& config, const DataSim& In);
  ~FullRecoTask();
  
  //int EventLoop(THypHi_Event *event,std::vector<TUTracker_Event*> *UTrackerEvents,Ana_Event* OutTree);

#ifdef ROOT6
  int EventLoop(const TG4Sol_Event& ev, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, MCAnaEventG4Sol* OutTree);
#else
  int EventLoop(const TG4Sol_Event& ev, const std::vector<TClonesArray*>& hits, MCAnaEventG4Sol* OutTree);
#endif  
  void AttachHisto(Ana_Hist* h);

private :

  const THyphiAttributes Attributes;

  FullRecoEvent REvent;
  Ana_Hist* AnaHisto;

  TDataBuilder* det_build;

  std::list<TDataProcess<FullRecoEvent,MCAnaEventG4Sol>*> list_processMC;

  int EventProcess(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree);
  

};

#endif
