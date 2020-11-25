#include "FullRecoTask.hh"

#include "TBuildDetectorLayerPlaneDAF.h"


#include "TBayesFinder.h"
//#include "TFinderCM.h"
#include "TKalmanFilter_DAF.h"
#include "CheckField.h"
#include "TCheckRZ.h"
//#include "TKalmanFilter_FRS.h"


/****************************************************************************************************/
/****************************************************************************************************/

// FullRecoTask::FullRecoTask():Attributes(),REvent()
// {
//   Attributes._logger->info(" *** > default FullRecoTask instance created !");
//   det_build = 0;
//   AnaHisto = 0;
//   list_processMC.resize(0);
// }



FullRecoTask::FullRecoTask(const FullRecoConfig& config, const DataSim& In):Attributes(config,In),REvent()
{
  Attributes._logger->info(" *** > FullRecoTask instance created | Reconstruction requested are :");
  det_build = nullptr;
  AnaHisto = nullptr;

  //det_build = new TBuildDetectorLayerPlane(Attributes,Attributes.beam_only);
  det_build = new TBuildDetectorLayerPlaneDAF(Attributes);
  //det_build = new TTestUnits(Attributes,"layerDAF");

  //list_process.push_back(new TKalmanFilter_DAF(Attributes) );
  list_processMC.emplace_back(new CheckField(Attributes));
  //list_processMC.emplace_back(new TBayesFinder(Attributes));
  list_processMC.emplace_back(new TCheckRZ(Attributes));
  //list_processMC.emplace_back(new TFinderCM(Attributes));
  list_processMC.emplace_back(new TKalmanFilter_DAF(Attributes) );

}


FullRecoTask::~FullRecoTask()
{
  //delete REvent;
  //FieldMap.clear();
  AnaHisto=nullptr;
  delete det_build;
  det_build = nullptr;
  for(auto process = list_processMC.begin(), process_last = list_processMC.end(); process!=process_last;++process)
    {
      delete (*process);
      *process = nullptr;
    }
}

void FullRecoTask::AttachHisto(Ana_Hist* h)
{
  AnaHisto = h;
  det_build->Init(h);
  for(auto process = list_processMC.begin(),process_last = list_processMC.end(); process!=process_last;++process)
    (*process)->Init(h);

}

//int FullRecoTask::EventLoop(THypHi_Event *event,std::vector<TUTracker_Event*> *UTrackerEvents,Ana_Event* OutTree)
#ifdef ROOT6
int FullRecoTask::EventLoop(const TG4Sol_Event& ev, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, MCAnaEventG4Sol* OutTree)
#else
int FullRecoTask::EventLoop(const TG4Sol_Event& ev, const std::vector<TClonesArray*>& hits, MCAnaEventG4Sol* OutTree)
#endif
{

  REvent.Clear();
  int return_build = (*dynamic_cast<TBuildDetectorLayerPlaneDAF*>(det_build))(ev,hits,REvent,OutTree);
  
  if(return_build !=0)
    return -1;

  int result_process = EventProcess(REvent,OutTree);

  return result_process;

}



int FullRecoTask::EventProcess(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree)
{

  for(auto process = list_processMC.begin(), process_last = list_processMC.end(); process!=process_last;++process)
    {
      //std::std::cout<<"courrent process:"<<(*process)->signature<<std::endl;
      if( (*(*process))(RecoEvent,OutTree) != 0)
	{
	  AnaHisto->h_task_exit.h->Fill((*process)->signature.c_str(),1);
	  Attributes._logger->warn("courrent process:{} failed", (*process)->signature);
	  return -1;
	}
    }
  return 0;
}







