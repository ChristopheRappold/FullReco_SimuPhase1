#include "FullRecoTask.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"
#include "TBuildDetectorLayerPlaneDAF.h"
#include "TBuildWASACalibrationLayerPlane.h"
#include "TBuildRestarter.h"

#include "TBayesFinder.h"
//#include "TFinderCM.h"
#include "TFindingPerf.h"
#include "TKalmanFilter_DAF.h"
#include "CheckField.h"
#include "TCheckFiberXUV.h"
#include "TCheckFiberTrack.h"
#include "TFragmentFinder.h"
#include "TCheckRZ.h"
#include "TFlatMCOutputML.h"
#include "TPrimaryVertex.h"
//#include "TPrimaryVertex_Si.h"
#include "TDecayVertex.h"
#include "TRiemannFinder.h"
#include "TFindingPerf.h"
#include <list>
//#include "TKalmanFilter_FRS.h"

/****************************************************************************************************/
/****************************************************************************************************/

//template<class TEOut> FullRecoTask<TEOut>::FullRecoTask():Attributes(),REvent()
// {
//   Attributes._logger->info(" *** > default FullRecoTask instance created !");
//   det_build = 0;
//   AnaHisto = 0;
//   list_processMC.resize(0);
// }



template<class TEOut>
FullRecoTask<TEOut>::FullRecoTask(const FullRecoConfig& config, const DataSimExp& In): Attributes(config,In),REvent()
{
  Attributes._logger->info(" *** > FullRecoTask instance created | Reconstruction requested are :");

  det_build = nullptr;
  AnaHisto = nullptr;

  if(Attributes.TaskConfig.Task_ReStart)
    det_build = new TBuildRestarter<TEOut>(Attributes);
  else
    {
      if constexpr(recotask::HasMC_Particle<TEOut>::value == true)
        det_build = new TBuildDetectorLayerPlaneDAF(Attributes);
      else
        det_build = new TBuildWASACalibrationLayerPlane(Attributes);
    }

  //list_process.push_back(new TKalmanFilter_DAF(Attributes) );
  // if(Attributes.TaskConfig.Task_CheckField)
  //   list_processMC.emplace_back(new CheckField<TEOut>(Attributes));
  // if(Attributes.TaskConfig.Task_PrimaryVtx)
  //   list_processMC.emplace_back(new TPrimaryVertex<TEOut>(Attributes));
  // if(Attributes.TaskConfig.Task_BayesFinder)
  //   list_processMC.emplace_back(new TBayesFinder<TEOut>(Attributes));
  // if(Attributes.TaskConfig.Task_RiemannFinder)
  //   list_processMC.emplace_back(new TRiemannFinder<TEOut>(Attributes));
  // //if(Attributes.TaskConfig.Task_FinderCM)
  // //  list_processMC.emplace_back(new TFinderCM(Attributes));
  // if(Attributes.TaskConfig.Task_FindingPerf)
  //   list_processMC.emplace_back(new TFindingPerf<TEOut>(Attributes));
  // if(Attributes.TaskConfig.Task_CheckRZ)
  //   list_processMC.emplace_back(new TCheckRZ<TEOut>(Attributes));
  // if(Attributes.TaskConfig.Task_KalmanDAF)
  //   list_processMC.emplace_back(new TKalmanFilter_DAF<TEOut>(Attributes));
  // if(Attributes.TaskConfig.Task_DecayVtx)
  //   list_processMC.emplace_back(new TDecayVertex<TEOut>(Attributes));
  // if(Attributes.TaskConfig.Task_FlatMCOutputML)
  //   list_processMC.emplace_back(new TFlatMCOutputML<TEOut>(Attributes));


  for(const auto& Tid : Attributes.TaskConfig.Task_Order)
    {
      switch(Tid)
        {
          case Task::TASKCHECKFIELD:
            if(Attributes.TaskConfig.Task_CheckField)
              list_processMC.emplace_back(new CheckField<TEOut>(Attributes));
            break;
          case Task::TASKPRIMARYVTX:
            if(Attributes.TaskConfig.Task_PrimaryVtx)
              list_processMC.emplace_back(new TPrimaryVertex<TEOut>(Attributes));
            break;
          //  case Task::TASKPRIMARYVTX_SI:
          //    if(Attributes.TaskConfig.Task_PrimaryVtx_Si)
          //      list_processMC.emplace_back(new TPrimaryVertex_Si<TEOut>(Attributes));
          //    break;
          case Task::TASKFLATMCOUTPUTML:
            if(Attributes.TaskConfig.Task_FlatMCOutputML)
              list_processMC.emplace_back(new TFlatMCOutputML<TEOut>(Attributes));
            break;
          case Task::TASKCHECKFIBERXUV:
            if(Attributes.TaskConfig.Task_CheckFiberXUV)
              list_processMC.emplace_back(new TCheckFiberXUV<TEOut>(Attributes));
            break;
          case Task::TASKCHECKFIBERTRACK:
            if(Attributes.TaskConfig.Task_CheckFiberTrack)
              list_processMC.emplace_back(new TCheckFiberTrack<TEOut>(Attributes));
            break;
          case Task::TASKFRAGMENTFINDER:
            if(Attributes.TaskConfig.Task_FragmentFinder)
              list_processMC.emplace_back(new TFragmentFinder<TEOut>(Attributes));
            break;
          case Task::TASKBAYESFINDER:
            if(Attributes.TaskConfig.Task_BayesFinder)
              list_processMC.emplace_back(new TBayesFinder<TEOut>(Attributes));
            break;
          case Task::TASKRIEMANNFINDER:
            if(Attributes.TaskConfig.Task_RiemannFinder)
              list_processMC.emplace_back(new TRiemannFinder<TEOut>(Attributes));
            break;
          case Task::TASKFINDERCM:
            //if(Attributes.TaskConfig.Task_FinderCM)
            //  list_processMC.emplace_back(new TFinderCM(Attributes));
            break;
          case Task::TASKFINDINGPERF:
            if(Attributes.TaskConfig.Task_FindingPerf)
              list_processMC.emplace_back(new TFindingPerf<TEOut>(Attributes));
            break;
          case Task::TASKCHECKRZ:
            if(Attributes.TaskConfig.Task_CheckRZ)
              list_processMC.emplace_back(new TCheckRZ<TEOut>(Attributes));
          case Task::TASKKALMANDAF:
            if(Attributes.TaskConfig.Task_KalmanDAF)
              list_processMC.emplace_back(new TKalmanFilter_DAF<TEOut>(Attributes));
            break;
          case Task::TASKDECAYVTX:
            if(Attributes.TaskConfig.Task_DecayVtx)
              list_processMC.emplace_back(new TDecayVertex<TEOut>(Attributes));
            break;
          default:
            break;
        }
    }

  for(auto task : list_processMC)
    Attributes._logger->info(" -> Task : {}",task->signature);
}


template<class TEOut>
FullRecoTask<TEOut>::~FullRecoTask()
{
  AnaHisto=nullptr;
  delete det_build;
  det_build = nullptr;
  for(auto process = list_processMC.begin(), process_last = list_processMC.end(); process!=process_last;++process)
    {
      delete (*process);
      *process = nullptr;
    }

  Attributes.SaveToDatabase();
}

template<class TEOut>
void FullRecoTask<TEOut>::AttachHisto(Ana_Hist* h)
{
  AnaHisto = h;
  det_build->Init(h);
  for(auto process = list_processMC.begin(),process_last = list_processMC.end(); process!=process_last;++process)
    (*process)->Init(h);
}

template<class TEOut>
void FullRecoTask<TEOut>::SetEventMetadata(AnaEvent_Metadata& metadata)
{
  metadata.NameIn = Attributes.NameIn;
  metadata.NameOut = Attributes.NameOut;
  metadata.DateOfRun = Attributes.DateOfRun;
  metadata.Hash = Attributes.Hash;
  metadata.FirstStep = det_build->signature;
  metadata.FinalStep = list_processMC.size() == 0 ? "None" : list_processMC.back()->signature;
  metadata.G4_simu = Attributes.G4_simu;
  metadata.NEvent = Attributes.NEvent;
  metadata.StartEvent = Attributes.Config.Get<uint>("Start_Event");
  metadata.StopEvent = Attributes.Config.Get<uint>("Stop_Event");
  metadata.Nb_Fraction = Attributes.Nb_Fraction;
  metadata.Wasa_Side = Attributes.Wasa_Side;
  metadata.Wasa_FieldMap = Attributes.Wasa_FieldMap;
  metadata.Field_Strength = Attributes.Field_Strength;
  metadata.Wasa_FieldMapName = Attributes.Wasa_FieldMapName;

  // Change ? -> Include event by event MetaData or full file metadata

}

//int template<class TEOut> FullRecoTask<TEOut>::EventLoop(THypHi_Event *event,std::vector<TUTracker_Event*> *UTrackerEvents,Ana_Event* OutTree)
#ifdef ROOT6
template<class TEOut>
int FullRecoTask<TEOut>::EventLoop(const TG4Sol_Event& ev, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, TEOut* OutTree)
#else
template<class TEOut>
int FullRecoTask<TEOut>::EventLoop(const TG4Sol_Event& ev, const std::vector<TClonesArray*>& hits, TEOut* OutTree)
#endif
{

  REvent.Clear();
  int return_build = (*dynamic_cast<TBuildDetectorLayerPlaneDAF*>(det_build))(ev,hits,REvent,OutTree);
  
  if(return_build !=0)
    return -1;

  int result_process = EventProcess(REvent,OutTree);

  return result_process;
}

#ifdef ROOT6
template<class TEOut>
int FullRecoTask<TEOut>::EventLoop(const MCAnaEventG4Sol& RestartEvent, TEOut* OutTree)
#else
template<class TEOut>
int FullRecoTask<TEOut>::EventLoop(TEOut* RestartEvent, TEOut* OutTree)
#endif
{

  REvent.Clear();
  int return_build = (*dynamic_cast<TBuildRestarter<TEOut>*>(det_build))(RestartEvent,REvent,OutTree);

  if(return_build !=0)
    return -1;

  int result_process = EventProcess(REvent,OutTree);

  return result_process;
}

#ifdef ROOT6
template<class TEOut>
int FullRecoTask<TEOut>::EventLoop(const EventWASAUnpack& UnpackEvent, TEOut* OutTree)
#else
template<class TEOut>
int FullRecoTask<TEOut>::EventLoop(const EventWASAUnpack& UnpackEvent, TEOut* OutTree)
#endif
{

  REvent.Clear();
  int return_build = (*dynamic_cast<TBuildWASACalibrationLayerPlane*>(det_build))(UnpackEvent,REvent,OutTree);

  if(return_build !=0)
    return -1;

  int result_process = EventProcess(REvent,OutTree);

  return result_process;
}


template<class TEOut>
int FullRecoTask<TEOut>::EventProcess(FullRecoEvent& RecoEvent,TEOut* OutTree)
{

  for(auto process = list_processMC.begin(), process_last = list_processMC.end(); process!=process_last;++process)
    {
      //std::std::cout<<"courrent process:"<<(*process)->signature<<std::endl;
      if( (*(*process))(RecoEvent,OutTree) != 0)
	{
	  AnaHisto->h_task_exit.h->Fill((*process)->signature.c_str(),1);
	  Attributes._logger->warn("current process:{} failed", (*process)->signature);
	  return -1;
	}
    }
  return 0;
}

template class FullRecoTask<MCAnaEventG4Sol>;
template class FullRecoTask<Ana_WasaEvent>;
