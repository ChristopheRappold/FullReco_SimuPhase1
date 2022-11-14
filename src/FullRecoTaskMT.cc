#include "FullRecoTaskMT.hh"

//#include "CheckField.h"
#include "TBuildDetectorLayerPlaneDAF_MT.h"
#include "TKalmanFilter_DAF_MT.h"
#include "TMergerOutput_MT.h"
//#include "TKalmanFilter_FRS.h"
#include "MTQueue.hh"
#include <thread>

/****************************************************************************************************/
/****************************************************************************************************/

// FullRecoTaskMT::FullRecoTaskMT():Attributes(),REvent()
// {
//   Attributes._logger->info(" *** > default FullRecoTaskMT instance created !");
//   det_build = 0;
//   AnaHisto = 0;
//   list_processMC.resize(0);
// }

FullRecoTaskMT::FullRecoTaskMT(const FullRecoConfig& conf, const DataSimExp& In)
  : Attributes(conf, In), REvent(Attributes.MTsetting.NQueue)
{
  Attributes._logger->info(" *** > FullRecoTaskMT instance created | Reconstruction requested are :");
  //det_build = nullptr;
  AnaHisto  = nullptr;
  //merger_out = nullptr;

  // for(size_t i = 0; i < Attributes.NthreadsQueue; ++i)
  //   REvent.emplace_back(FullRecoEvent(i));

  // REvent.resize(Attributes.NthreadsQueue);

  // det_build = new TBuildDetectorLayerPlane(Attributes,Attributes.beam_only);
  for(int i=0;i<Attributes.MTsetting.NBuilder;++i)
    list_det_build.emplace_back(new TBuildDetectorLayerPlaneDAF_MT(Attributes));
  // det_build = new TTestUnits(Attributes,"layerDAF");

  // list_process.push_back(new TKalmanFilter_DAF(Attributes) );
  //list_processMC.emplace_back(new CheckField(Attributes));
  //list_processMC.emplace_back(new TKalmanFilter_DAF(Attributes));
  for(int i =0 ;i<Attributes.MTsetting.NKalman;++i)
    list_processMC_MT.emplace_back(new TKalmanFilter_DAF_MT(Attributes));

  for(int i =0 ;i<Attributes.MTsetting.NMerger;++i)
    list_merger_out.emplace_back( new TMergerOutput_MT(Attributes));
}

FullRecoTaskMT::~FullRecoTaskMT()
{
  // delete REvent;
  // FieldMap.clear();
  AnaHisto = nullptr;
  
  for(auto process = list_det_build.begin(), process_last = list_det_build.end(); process != process_last;      ++process)
    {
      delete (*process);
      *process = nullptr;
    }
  for(auto process = list_processMC_MT.begin(), process_last = list_processMC_MT.end(); process != process_last;
      ++process)
    {
      delete(*process);
      *process = nullptr;
    }
  for(auto process = list_merger_out.begin(), process_last = list_merger_out.end(); process != process_last;      ++process)
    {
      delete (*process);
      *process = nullptr;
    }
}

void FullRecoTaskMT::AttachHisto(Ana_Hist* h)
{
  AnaHisto = h;
  for(auto& process : list_det_build)
    process->Init(h);

  for(auto process = list_processMC_MT.begin(), process_last = list_processMC_MT.end(); process != process_last;
      ++process)
    {
      (*process)->Init(h);
    }
  for(auto& process : list_merger_out)
    process->Init(h);
  
}

int FullRecoTaskMT::Run(Long64_t startEvent, Long64_t stopEvent, TTree* InTree, const TG4Sol_Event& ev,
                        const std::vector<TClonesArray*>& hits, TTree* OutTree, MCAnaEventG4Sol* OutEvent)
{
  using id_th = int;
  Attributes._logger->debug("!> MT Run start: {}",Attributes.MTsetting.NQueue);

  fast_blocking_queue<std::tuple<id_th, ReturnRes::InfoM> > queueInit(Attributes.MTsetting.NQueue);
  for(int i = 0; i < Attributes.MTsetting.NQueue; ++i)
    {
      Attributes._logger->debug("!> queueInit : {}",i); 
      REvent[i].idThread = i;
      queueInit.push(std::make_tuple(i, ReturnRes::Fine));
    }
  Attributes._logger->debug("!> queueInit size :{}",queueInit.empty());
  for(size_t i = 0; i<REvent.size();++i)
    Attributes._logger->debug("!> RecoEvent#{} -> id {} ",i,REvent[i].idThread);
  
  fast_blocking_queue<std::tuple<id_th, ReturnRes::InfoM> > queueStep1(Attributes.MTsetting.NQueue);
  fast_blocking_queue<std::tuple<id_th, ReturnRes::InfoM> > queueStep2(Attributes.MTsetting.NQueue);
  fast_blocking_queue<std::tuple<id_th, ReturnRes::InfoM> > queueStep3(Attributes.MTsetting.NQueue);

  const int sizeStep1 = Attributes.MTsetting.NBuilder;
  const int sizeStep2 = Attributes.MTsetting.NKalman;
  const int sizeStep3 = Attributes.MTsetting.NMerger;

  struct nextingEvent {
    const Long64_t Start;
    const Long64_t Stop;
    std::atomic_ullong Current;
    nextingEvent(Long64_t f, Long64_t l, Long64_t c):Start(f),Stop(l) { Current = c;}  
  };
  nextingEvent nextEvent(startEvent, stopEvent, startEvent);
  auto& logging = Attributes._logger;
 
  auto f_nextEvent = [&](nextingEvent& nEvent, std::vector<FullRecoEvent>& REvent) {

    logging->debug("MT0> start thread nextEvent / queueInit: {}",queueInit.empty());

    while(1)
      {
	std::tuple<id_th, ReturnRes::InfoM> Status;
	queueInit.pop(Status);
	auto currentEvent = nEvent.Current.fetch_add(1, std::memory_order_relaxed);
	logging->debug("MT0> Current Event:{}",currentEvent);
	if(currentEvent < nEvent.Stop)
	  {
	    assert(std::get<0>(Status) == REvent[std::get<0>(Status)].idThread);
	    REvent[std::get<0>(Status)].idEvent = currentEvent;
	    queueStep1.push(Status);
	  }
	else
	  {
	    logging->debug("MT0> finishing !  init {} / q1 : {} / q2 : {} / q3 : {}",queueInit.empty(), queueStep1.empty(), queueStep2.empty(), queueStep3.empty() );
	    std::get<0>(Status) = sizeStep1;
	    std::get<1>(Status) = ReturnRes::EndRun;
	    queueStep1.push(Status);
	    return; 
	  }
      }
  };
  
  auto& list_builder = list_det_build;

  auto f_det_build = [&](int idTh, TTree* InTree, const TG4Sol_Event& ev,
                         const std::vector<TClonesArray*>& hits, std::vector<FullRecoEvent>& REvent) {

    logging->debug("MT1> start thread det_build ");

    while(1)
      {
	std::tuple<id_th, ReturnRes::InfoM> Status;
        queueStep1.pop(Status);
	if(std::get<1>(Status) == ReturnRes::EndRun)
	  {
	    int tempS = std::get<0>(Status);
	    logging->debug("MT1> finishing ! | init {} Q1 {} Q2 {} Q3 {} : {}",queueInit.empty(), queueStep1.empty(), queueStep2.empty(), queueStep3.empty(), tempS );
	    if(tempS == sizeStep1)
	      {
		auto Status2 = Status;
		std::get<0>(Status2) = sizeStep2;
		queueStep2.push(Status2);
		logging->debug("MT1> finishing ! | push Q2 finished status {}",tempS);
	      }
	    --tempS;
	    std::get<0>(Status) = tempS;
	    if(tempS != 0)
	      queueStep1.push(Status);
	    return;
	  }
	
        assert(std::get<0>(Status) == REvent[std::get<0>(Status)].idThread);

        auto idEvent = REvent[std::get<0>(Status)].idEvent;

	InTree->GetEntry(idEvent);
        ReturnRes::InfoM return_build =
	  (*dynamic_cast<TBuildDetectorLayerPlaneDAF_MT*>(list_builder[idTh]))(ev, hits, std::ref(REvent[std::get<0>(Status)]));

	std::get<1>(Status) = return_build;
        queueStep2.push(Status);
      }
  };

  auto& list_processMC_RK = list_processMC_MT;

  auto f_RK_daf = [&](int htId, std::vector<FullRecoEvent>& RecoEvent) {

    std::stringstream ss1;
    ss1 <<  std::this_thread::get_id() ;
    logging->debug("MT2> start thread RK#{} {}",htId,ss1.str());
    list_processMC_RK[htId]->InitMT();
    while(1)
      {
        std::tuple<id_th, ReturnRes::InfoM> Status;
        queueStep2.pop(Status);
	if(std::get<1>(Status) == ReturnRes::EndRun)
	  {
	    int tempS = std::get<0>(Status);
	    logging->debug("MT2> finishing ! | init {} Q1 {} Q2 {} Q3 {} : {}",queueInit.empty(), queueStep1.empty(), queueStep2.empty(), queueStep3.empty(), tempS );
	    if(tempS == sizeStep1)
	      {
		auto Status2 = Status;
		std::get<0>(Status2) = sizeStep3;
		queueStep3.push(Status2);
		logging->debug("MT2> finishing ! | push Q3 finished status {}",tempS);
	      }
	    --tempS;
	    std::get<0>(Status) = tempS;
	    if(tempS != 0)
	      queueStep2.push(Status);

	    return;
	  }
	else if(std::get<1>(Status) != ReturnRes::Fine)
          queueStep3.push(Status);

        assert(std::get<0>(Status) == REvent[std::get<0>(Status)].idThread);

        ReturnRes::InfoM return_process = (*list_processMC_RK[htId])(std::ref(RecoEvent[std::get<0>(Status)]),nullptr);

	std::get<1>(Status)  = return_process;
        queueStep3.push(Status);
      }
  };

  auto& list_merger = list_merger_out;
  
  auto f_Output = [&](int idTh, std::vector<FullRecoEvent>& RecoEvent, TTree* OutTree, MCAnaEventG4Sol* OutEvent) {
    logging->debug("MT3> start thread Output");

    while(1)
      {

        std::tuple<id_th, ReturnRes::InfoM> Status;
        queueStep3.pop(Status);
	if(std::get<1>(Status) == ReturnRes::EndRun)
	  {
	    int tempS = std::get<0>(Status);
	    logging->debug("MT3> finishing ! | init {} Q0 {} Q1 {} Q3 {} : {}",queueInit.empty(), queueStep1.empty(), queueStep2.empty(), queueStep3.empty(), tempS );
	    --tempS;
	    std::get<0>(Status) = tempS;
	    if(tempS != 0)
	      queueStep3.push(Status);
	    return;
	  }
        assert(std::get<0>(Status) == REvent[std::get<0>(Status)].idThread);

        (*dynamic_cast<TMergerOutput_MT*>(list_merger[idTh]))(std::ref(RecoEvent[std::get<0>(Status)]), std::get<1>(Status), OutEvent);
	
	OutTree->Fill();
	OutEvent->Clear();
	
        RecoEvent[std::get<0>(Status)].Clear();
        std::get<1>(Status) = ReturnRes::Fine;
        queueInit.push(Status);

      }
  };

  std::thread t0_nE(f_nextEvent,std::ref(nextEvent),std::ref(REvent));
  std::vector<std::thread> t1s;
  for(int i=0;i<sizeStep1;++i)
    t1s.emplace_back(f_det_build,i,std::ref(InTree),std::ref(ev),std::ref(hits),std::ref(REvent));

  std::vector<std::thread> t2s;
  for(int i=0;i<sizeStep2;++i)
    t2s.emplace_back(f_RK_daf,i,std::ref(REvent));
  // std::thread t3_RK1(f_RK_daf,1,std::ref(REvent));
  //std::thread t3_RK2(f_RK_daf,2,std::ref(REvent));
  std::vector<std::thread> t3s;
  for(int i=0;i<sizeStep3;++i)
    t3s.emplace_back(f_Output,i,std::ref(REvent),OutTree,OutEvent);

  // t1_nE.join();
  // t2_build.join();
  // t3_RK0.join();
  // t3_RK1.join();
  // t3_RK2.join();
  // t4_out.join();

  t0_nE.join();
  for(auto& t1 : t1s)
    t1.join();
  for(auto& t2 : t2s)
    t2.join();
  for(auto& t3 : t3s)
    t3.join();
  
  return 0;
}
