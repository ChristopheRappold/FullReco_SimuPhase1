#include "FullRecoTaskZMQ.hh"

//#include "CheckField.h"
// #include "TBuildDetectorLayerPlaneDAF_ZMQ.h"
// #include "TKalmanFilter_DAF_ZMQ.h"
// #include "TMergerOutput_ZMQ.h"
//#include "TKalmanFilter_FRS.h"
//#include "MTQueue.hh"
#include "ThreadTasks.hh"
#include <thread>

using namespace ZMQ;

/****************************************************************************************************/
/****************************************************************************************************/

// FullRecoTaskZMQ::FullRecoTaskZMQ():Attributes(),REvent()
// {
//   Attributes._logger->info(" *** > default FullRecoTaskZMQ instance created !");
//   det_build = 0;
//   AnaHisto = 0;
//   list_processMC.resize(0);
// }

FullRecoTaskZMQ::FullRecoTaskZMQ(const FullRecoConfig& conf, const DataSim& In) : Attributes(conf, In)
{
  Attributes._logger->info(" *** > FullRecoTaskZMQ instance created | Reconstruction requested are :");
  // det_build = nullptr;
  AnaHisto = nullptr;
  // merger_out = nullptr;

  // for(size_t i = 0; i < Attributes.NthreadsQueue; ++i)
  //   REvent.emplace_back(FullRecoEvent(i));

  // REvent.resize(Attributes.NthreadsQueue);

  // det_build = new TBuildDetectorLayerPlane(Attributes,Attributes.beam_only);
  for(int i = 0; i < Attributes.NBuilder; ++i)
    list_det_build.emplace_back(new TBuildDetectorLayerPlaneDAF_ZMQ(Attributes));
  // det_build = new TTestUnits(Attributes,"layerDAF");

  // list_process.push_back(new TKalmanFilter_DAF(Attributes) );
  // list_processMC.emplace_back(new CheckField(Attributes));
  // list_processMC.emplace_back(new TKalmanFilter_DAF(Attributes));
  for(int i = 0; i < Attributes.NKalman; ++i)
    list_processMC_MT.emplace_back(new TKalmanFilter_DAF_ZMQ(Attributes));

  for(int i = 0; i < Attributes.NMerger; ++i)
    list_merger_out.emplace_back(new TMergerOutput_ZMQ(Attributes));

  IsMain      = Attributes.IsMain;
  sizeBuilder = Attributes.NBuilder;
  sizeKalman  = Attributes.NKalman;
  sizeMerger  = Attributes.NMerger;

  addr_initEvent = Attributes.addr_initEvent;

  addr_frontBuilder = Attributes.addr_frontBuilder;
  addr_backFitter   = Attributes.addr_backFitter;

  addr_frontFitter = Attributes.addr_frontFitter;
  addr_backMerger  = Attributes.addr_backMerger;

  addr_frontMerger = Attributes.addr_frontMerger;
  addr_backEnd     = Attributes.addr_backEnd;

  addr_control = Attributes.addr_control;

  addr_monitor = Attributes.addr_monitor;

  Attributes._logger->info("IsMain :{} / B: {} / K: {} / M: {}", IsMain, sizeBuilder, sizeKalman, sizeMerger);
  Attributes._logger->info("Addresses :");
  Attributes._logger->info("nEvent   : {}", addr_initEvent);
  Attributes._logger->info("F_Builder: {}",addr_frontBuilder);
  Attributes._logger->info("B_Fitter : {}",addr_backFitter);
  Attributes._logger->info("F_Fitter : {}",addr_frontFitter);
  Attributes._logger->info("B_Merger : {}",addr_backMerger);
  Attributes._logger->info("F_Merger : {}",addr_frontMerger);
  Attributes._logger->info("End      : {}",addr_backEnd);
  Attributes._logger->info("control  : {}",addr_control);
  Attributes._logger->info("monitor  : {}",addr_monitor);
}

FullRecoTaskZMQ::~FullRecoTaskZMQ()
{
  // delete REvent;
  // FieldMap.clear();
  AnaHisto = nullptr;

  for(auto process = list_det_build.begin(), process_last = list_det_build.end(); process != process_last; ++process)
    {
      delete(*process);
      *process = nullptr;
    }
  for(auto process = list_processMC_MT.begin(), process_last = list_processMC_MT.end(); process != process_last;
      ++process)
    {
      delete(*process);
      *process = nullptr;
    }
  for(auto process = list_merger_out.begin(), process_last = list_merger_out.end(); process != process_last; ++process)
    {
      delete(*process);
      *process = nullptr;
    }
}

void FullRecoTaskZMQ::AttachHisto(Ana_Hist* h)
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

int FullRecoTaskZMQ::Run(Long64_t startEvent, Long64_t stopEvent, TTree* InTree, const TG4Sol_Event& ev,
                         const std::vector<TClonesArray*>& hits, TTree* OutTree, MCAnaEventG4Sol* OutEvent)
{
  Attributes._logger->debug("!> MT Run start: {}", Attributes.NQueue);

  
  std::thread t0_nE;
  std::thread t1_proxy;
  std::thread t2_proxy;
  std::thread t3_proxy;
  std::thread t4_end;

  nextingEvent nextEvent(startEvent, stopEvent, startEvent,20);

  
  if(IsMain)
    {

      c_nextEvent f_nextEvent(contextTask,Attributes._logger);
  
      c_proxy f_proxy(contextTask);
      
      // auto f_proxy = [&context](const std::string& addr_F, const std::string& addr_B, const std::string& addr_C,
      //                           const std::string& nameProxy) {
      //   proxyPullPushNN proxy(context, std::ref(addr_F), std::ref(addr_B), std::ref(addr_C), "Console",
      //                         std::ref(nameProxy));
      //   proxy.Run();
      // };
      c_finalEvent<1000,100> f_finalEvent(contextTask,Attributes._logger);

      t0_nE = std::thread(f_nextEvent, std::ref(nextEvent), addr_initEvent, addr_control, addr_monitor);

      t1_proxy = std::thread(f_proxy, addr_frontBuilder, addr_backFitter, addr_control, "BuilderFitter");
      t2_proxy = std::thread(f_proxy, addr_frontFitter, addr_backMerger, addr_control, "FitterMerger");
      t3_proxy = std::thread(f_proxy, addr_frontMerger, addr_backEnd, addr_control, "MergerEnd");

      t4_end = std::thread(f_finalEvent, stopEvent - startEvent, addr_backEnd, addr_control, addr_initEvent, addr_monitor);


      // std::thread t0_nE(f_nextEvent,std::ref(nextEvent));

  // std::thread t1_proxy(f_proxyBuilderFitter);
  // std::thread t2_proxy(f_proxyFitterMerger);
  // std::thread t3_proxy(f_proxyMergerEnd);

  // std::thread t4_end(f_finalEvent, stopEvent-startEvent);
    }

  std::vector<std::thread> t1s;

  if(sizeBuilder>0)
    {
      c_det_build f_det_build(contextTask, Attributes._logger, list_det_build);
      
      for(int i = 0; i < sizeBuilder; ++i)
	t1s.emplace_back(f_det_build, i, std::ref(InTree), std::ref(ev), std::ref(hits), addr_initEvent, addr_frontBuilder,
			 addr_control, addr_monitor);

    }

  std::vector<std::thread> t2s;

   if(sizeKalman>0)
     {
       c_RK f_RK(contextTask, Attributes._logger, list_processMC_MT);
       
       for(int i = 0; i < sizeKalman; ++i)
	 t2s.emplace_back(f_RK, i, addr_backFitter, addr_frontFitter, addr_control, addr_monitor);

     }

  std::vector<std::thread> t3s;

  if(sizeMerger>0)
    {
      c_Output f_Output(contextTask, Attributes._logger, list_merger_out);
  
      for(int i = 0; i < sizeMerger; ++i)
	t3s.emplace_back(f_Output, i, OutTree, OutEvent, addr_backMerger, addr_frontMerger, addr_control, addr_monitor);

    }
  // std::thread t0_nE(f_nextEvent,std::ref(nextEvent));

  // std::thread t1_proxy(f_proxyBuilderFitter);
  // std::thread t2_proxy(f_proxyFitterMerger);
  // std::thread t3_proxy(f_proxyMergerEnd);

  // std::thread t4_end(f_finalEvent, stopEvent-startEvent);

  // t1_nE.join();
  // t2_build.join();
  // t3_RK0.join();
  // t3_RK1.join();
  // t3_RK2.join();
  // t4_out.join();
  if(IsMain)
    {
      t0_nE.join();

      t1_proxy.detach();
      t2_proxy.detach();
      t3_proxy.detach();
    }
  for(auto& t1 : t1s)
    t1.join();

  // std::this_thread::sleep_for(std::chrono::milliseconds(200));
  for(auto& t2 : t2s)
    t2.join();

  // std::this_thread::sleep_for(std::chrono::milliseconds(200));

  for(auto& t3 : t3s)
    t3.join();

  // std::this_thread::sleep_for(std::chrono::milliseconds(200));

  if(IsMain)
    t4_end.join();

  return 0;
}
