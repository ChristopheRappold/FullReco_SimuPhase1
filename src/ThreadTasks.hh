#ifndef THREADTASKS_H
#define THREADTASKS_H

#include "FullRecoEventZMQ.hh"
#include "TBuildDetectorLayerPlaneDAF_ZMQ.h"
#include "TKalmanFilter_DAF_ZMQ.h"
#include "TMergerOutput_ZMQ.h"
#include "TTree.h"
#include "ThreadModelZMQ.hh"

#include <chrono>
#include <thread>

using namespace ZMQ;

struct nextingEvent
{
  const Long64_t Start;
  const Long64_t Stop;
  const Long64_t Step;
  std::atomic_ullong Current;
  nextingEvent(Long64_t first, Long64_t last, Long64_t c, Long64_t step) : Start(first), Stop(last), Step(step)
  {
    Current = c;
  }
};

template <int Nstats = 100>
struct c_nextEvent
{

  zmq::context_t& context;
  const std::shared_ptr<spdlog::logger>& logging;

  explicit c_nextEvent(zmq::context_t& co, const std::shared_ptr<spdlog::logger>& l) : context(co), logging(l) {}

  void operator()(nextingEvent& nEvent, const std::string& addr_I, const std::string& addr_C, const std::string& addr_M)
  {

    zmq::socket_t sock_init(context, zmq::socket_type::push);
    sock_init.bind(addr_I.c_str()); //"tcp://127.0.0.1:10001");

    zmq::socket_t control(context, zmq::socket_type::sub);
    control.connect(addr_C.c_str());
    //control.setsockopt(ZMQ_SUBSCRIBE, "", 0);
    control.set(zmq::sockopt::subscribe, "");

    static const std::string tag_stats("stats_nEvent");
    zmq::socket_t monitor(context, zmq::socket_type::pub);
    monitor.connect(addr_M.c_str());

    std::this_thread::sleep_for(std::chrono::milliseconds(200));

    States state        = init;
    int count           = 0;
    auto time_ref_start = std::chrono::high_resolution_clock::now();
    auto pidSys = ::getpid(); //std::hash<std::thread::id>()( std::this_thread::get_id() );

    while(1)
      {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        EventStatus Status;

        switch(state)
          {
          case init:
            {
              logging->info("MT_init> start thread nextEvent");
              state = process;
            }
          case process:
            {
              // queue_init.pop(Status);
              auto currentEvent = nEvent.Current.fetch_add(nEvent.Step, std::memory_order_relaxed);

              logging->info("MTinit> Current Event:{}", currentEvent);
              if(currentEvent < nEvent.Stop)
                {
                  Status.status     = ReturnRes::Fine;
                  Status.BeginEvent = currentEvent;
                  Status.EndEvent = currentEvent + nEvent.Step > nEvent.Stop ? nEvent.Stop : currentEvent + nEvent.Step;
                  logging->info("MTinit> Current Range Event:[{}, {}]", Status.BeginEvent, Status.EndEvent);
                  try
                    {
                      msgpack::sbuffer sbuf;
                      msgpack::pack(sbuf, Status);
                      // zmq::message_t body_msg(sbuf.size());
                      zmq::const_buffer body_msg(sbuf.data(), sbuf.size());
                      // std::memcpy(body_msg.data(), sbuf.data(), sbuf.size());

                      sock_init.send(body_msg, zmq::send_flags::dontwait);
                    }
                  catch(zmq::error_t& e)
                    {
                      logging->error("Error catched in MTinit event<max ! {} {}", e.what(), e.num());
                    }
                  ++count;
                  if(count == Nstats)
                    {
                      state = stats;
                      count = 0;
                    }
                }
              else
                {
                  state = finishing;
                }
              break;
            }
          case stats:
            {
              auto time_stat = std::chrono::high_resolution_clock::now();
              std::chrono::duration<double> time_diff = time_stat - time_ref_start;

	      std::tuple<int,int,int,double> time_data = std::make_tuple(pidSys, 1, Nstats, time_diff.count());
              publish(monitor, tag_stats, time_data);

              state = process;
              break;
            }
          case finishing:
            {
              logging->info("MTinit> finishing !");
              //EventStatus Status;
              Status.status     = ReturnRes::EndRun;
              Status.BeginEvent = 0;
              Status.EndEvent   = 1;
              try
                {
                  msgpack::sbuffer sbuf;
                  msgpack::pack(sbuf, Status);

                  zmq::const_buffer body_msg(sbuf.data(), sbuf.size());

                  sock_init.send(body_msg,
                                 zmq::send_flags::dontwait); // zmq::buffer(Status), zmq::send_flags::dontwait);
                }
              catch(zmq::error_t& e)
                {
                  logging->error("Error catched in MTinit finishing ! {} {}", e.what(), e.num());
                }

              sock_init.unbind(addr_I.c_str());
              state = last;
              break;
            }
          case last:
            {
              logging->info("MTinit> last ! ");
              zmq::message_t body_msg;
              auto ret = control.recv(body_msg);
	      if(!ret)
		logging->warn("MTinit> last : no msg !");
	      
	      const std::string str_msg(static_cast<const char*>(body_msg.data()), body_msg.size());
              static const std::string terminate = "TERMINATE";

              if(str_msg == terminate)
                state = end;

              break;
            }
          case end:
            {
              logging->info("MTinit> End !");
              return;
            }
          }
      }
  }
};

struct c_proxy
{
  zmq::context_t& context;
  explicit c_proxy(zmq::context_t& c) : context(c) {}
  void operator()(const std::string& addr_F, const std::string& addr_B, const std::string& addr_C,
                  const std::string& nameProxy)
  {
    proxyPullPushNN proxy(context, std::ref(addr_F), std::ref(addr_B), std::ref(addr_C), "Console",
                          std::ref(nameProxy));
    proxy.Run();
  }
};

template <size_t N = 1000, int Nstats = 100>
struct c_finalEvent
{

  zmq::context_t& context;
  const std::shared_ptr<spdlog::logger>& logging;
  explicit c_finalEvent(zmq::context_t& c, const std::shared_ptr<spdlog::logger>& l) : context(c), logging(l) {}

  void operator()(int totalEvent, const std::string& addr_B, const std::string& addr_C, const std::string& addr_I,
                  const std::string& addr_M)
  {
    //Long64_t countEvent = 0;
    //Long64_t rangeEvent = 0;

    //const Long64_t baseCountEvent = N;

    std::array<bool, N> statusEvents;
    std::array<bool, N> NextstatusEvents;
    statusEvents.fill(false);
    NextstatusEvents.fill(false);

    std::set<Long64_t> missingEvents;

    // auto f_check = [](const std::array<bool, N>& s, int t) {
    //   std::set<Long64_t> missingEvents;
    //   for(size_t i = 0; i < t; ++i)
    //     {
    //       if(s.at(i) == false)
    //         missingEvents.insert(i);
    //     }
    //   return std::move(missingEvents);
    // };

    zmq::socket_t sockQ2_R(context, zmq::socket_type::pull);
    zmq::socket_t sockFinal(context, zmq::socket_type::pub);

    std::string addr_missing(addr_I);
    addr_missing += "9"; //"Missing";
    zmq::socket_t sockInit_back(context, zmq::socket_type::push);

    static const std::string tag_stats("stats_fEvent");
    zmq::socket_t monitor(context, zmq::socket_type::pub);

    try
      {
        // std::this_thread::sleep_for(std::chrono::milliseconds(100*idTh));

        sockQ2_R.connect(addr_B.c_str()); //("tcp://127.0.0.1:10001");
      }
    catch(zmq::error_t& e)
      {
        logging->error("bind Q2R ThQfinal {} : {}", 0, e.what());
        std::exit(-1);
      }
    try
      {
        // std::this_thread::sleep_for(std::chrono::milliseconds(100*idTh));
	monitor.connect(addr_M.c_str());
      }
    catch(zmq::error_t& e)
      {
        logging->error("bind monitor ThQfinal {} : {}", 0, e.what());
        std::exit(-1);
      }
    try
      {
        sockInit_back.bind(addr_missing.c_str());
      }
    catch(zmq::error_t& e)
      {
        logging->error("bind InitBack ThQfinal {} : {}", 0, e.what());
        std::exit(-1);
      }
    try
      {
        // std::this_thread::sleep_for(std::chrono::milliseconds(100*idTh));
        sockFinal.bind(addr_C.c_str());
      }
    catch(zmq::error_t& e)
      {
        logging->error("bind Final ThQfinal {} : {}", 0, e.what());
        std::exit(-1);
      }

    std::this_thread::sleep_for(std::chrono::milliseconds(200));

    States state        = init;
    int countStats      = 0;
    int previousStats = 0;
    auto time_ref_start = std::chrono::high_resolution_clock::now();
    auto pidSys = ::getpid(); //std::hash<std::thread::id>()( std::this_thread::get_id() );

    while(1)
      {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        switch(state)
          {
          case init:
            {
              logging->info("MT_final> start thread finalEvent");
              state = process;
            }
          case process:
            {
              EventStatus Status;

              zmq::message_t body_msg;
              try
                {
                  auto ret = sockQ2_R.recv(body_msg); //,zmq::recv_flags::none);//bufIn,zmq::recv_flags::dontwait);
		  if(!ret)
		    logging->warn("MT_final> process : no msg !");
		}
              catch(zmq::error_t& e)
                {
                  logging->error("Error catched in MTfinal#0 receiving ! {} {}", e.what(), e.num());
                }

              try
                {
                  msgpack::object_handle unpacked_body =
                      msgpack::unpack(static_cast<const char*>(body_msg.data()), body_msg.size());
                  msgpack::object obj = unpacked_body.get();
                  obj.convert(Status);
                }
              catch(msgpack::unpack_error& e)
                {
                  logging->error("msgpack :{}", e.what());
                  exit(-1);
                }
              // if(Status.status == ReturnRes::Fine)
              //   {

              //     Long64_t currentEvent  = Status.BeginEvent;
              //     Long64_t currentEventE = Status.EndEvent;
              //     if(currentEventE != currentEvent)
              //       logging->warn("MTfinal> received eventstatus as range not unique idEvent: B:{}, E:{}", currentEvent,
              //                     currentEventE);

              //     logging->info("MTfinal> Current Event:{}", currentEvent);
              //     if(currentEvent < rangeEvent * baseCountEvent)
              //       {
              //         auto it_miss = missingEvents.find(currentEvent);
              //         if(it_miss == missingEvents.end())
              //           logging->warn("MTfinal> missing Event {}: in passed ranged yet not missing event !",
              //                         currentEvent);
              //         else
              //           {
              //             logging->debug("MTfinal> missing Event {} found and removed !", currentEvent);
              //             missingEvents.erase(it_miss);
              //             if(missingEvents.empty())
              //               state = end;
              //           }
              //       }
              //     else if(currentEvent >= (rangeEvent + 1) * baseCountEvent)
              //       {
              //         std::set<Long64_t> temp_missing = f_check(statusEvents, baseCountEvent);

              //         for(auto& item : temp_missing)
              //           missingEvents.insert(item + rangeEvent * baseCountEvent);

              //         statusEvents.swap(NextstatusEvents);
              //         logging->debug("MTfinal> Higher range: {} [{}, {}] {} -> now {}", currentEvent, rangeEvent,
              //                        baseCountEvent, (rangeEvent + 1) * baseCountEvent,
              //                        currentEvent - (rangeEvent + 1) * baseCountEvent);
              //         currentEvent -= (rangeEvent + 1) * baseCountEvent;
              //         statusEvents[currentEvent] = true;
              //         ++countEvent;

              //         ++rangeEvent;
              //         NextstatusEvents.fill(false);
              //       }
              //     else
              //       {
              //         logging->debug("MTfinal> in range :{} -> now {}", currentEvent,
              //                        currentEvent - rangeEvent * baseCountEvent);

              //         auto it_miss = missingEvents.find(currentEvent);
              //         if(it_miss != missingEvents.end())
              //           {
              //             logging->debug("MTfinal> missing Event {} found and removed !", currentEvent);
              //             missingEvents.erase(it_miss);
              //             if(missingEvents.empty())
              //               state = end;
              //           }

              //         currentEvent -= rangeEvent * baseCountEvent;
              //         statusEvents[currentEvent] = true;
              //         ++countEvent;
              //       }
              //   }
              if(Status.status == ReturnRes::EndRun)
                {
                  state = finishing;
                }
              // else
              //   {
              //     int currentEvent  = Status.BeginEvent;
              //     int currentEvent2 = Status.EndEvent;
              //     logging->warn("MTfinal> Not fine / not endRun [Q{}]/ Current Event:[{}, {}]", Status.status,
              //                   currentEvent, currentEvent2);
              //   }

              ++countStats;
              if(countStats == Nstats)
                {
                  state      = stats;
                  countStats = 0;
                }

              break;
            }
          case stats:
            {
              auto time_stat = std::chrono::high_resolution_clock::now();
              std::chrono::duration<double> time_diff = time_stat - time_ref_start;
	      
	      std::tuple<int,int,int,double> time_data = std::make_tuple(pidSys, 1, Nstats*(previousStats+1), std::chrono::duration_cast<std::chrono::milliseconds>(time_diff).count());
              publish(monitor, tag_stats, time_data);
	      ++previousStats;
	      time_ref_start = time_stat;
              state = process;

              break;
            }
          case finishing:
            {
              logging->info("MTfinal> finishing !");
              state = last;
              break;
            }
          case end:
          case last:
            {
              logging->info("MTfinal> last ! "); // init {} / q0 : {} / q1 : {} / q2 : {}",queue_init.empty(),

              // std::set<Long64_t> last_missingEvents = f_check(statusEvents, totalEvent - rangeEvent * baseCountEvent);
              // for(auto& item : last_missingEvents)
              //   missingEvents.insert(item + rangeEvent * baseCountEvent);

              // std::stringstream ss1, ss2;
              // ss1 << "last: ";
              // for(size_t i = 0; i < baseCountEvent; ++i)
              //   ss1 << "[" << i << " " << statusEvents[i] << "] ";
              // logging->info(ss1.str());
              // logging->info("missing Event:");

              // for(auto e : missingEvents)
              //   ss2 << " " << e;

              // logging->info(ss2.str());

              if(missingEvents.empty())
                {
                  logging->info("MTfinal> sending TERMINATE !");
                  std::this_thread::sleep_for(std::chrono::milliseconds(10));
                  sockFinal.send(zmq::const_buffer("TERMINATE", 9));
                  std::this_thread::sleep_for(std::chrono::milliseconds(100));
                  return;
                }
              else
                {
                  logging->info("MTfinal> sending missing !");
                  logging->info("MTfinal> binding done to {}!", addr_missing);
                  for(auto idEvent : missingEvents)
                    {
                      try
                        {
                          EventStatus StatusO;
                          StatusO.status     = ReturnRes::Fine;
                          StatusO.BeginEvent = idEvent;
                          StatusO.EndEvent   = idEvent;

                          msgpack::sbuffer sbuf;
                          msgpack::pack(sbuf, StatusO);

                          zmq::const_buffer body_msg(sbuf.data(), sbuf.size());
                          logging->info("MTfinal> sending missing :{}", idEvent);
                          sockInit_back.send(body_msg, zmq::send_flags::dontwait);
                        }
                      catch(zmq::error_t& e)
                        {
                          logging->error("Error catched in MTfinal initback ! {} {}", e.what(), e.num());
                        }
                    }
                  logging->info("MTfinal> sending done ");

                  state = process;
                }
              break;
            }
          }
      }
  }
};

inline void f_builder(ZMQ::EventStatus& In, ZMQ::DataBuilderOut& Out, TBuildDetectorLayerPlaneDAF_ZMQ* process,
                      std::shared_ptr<spdlog::logger>& logging, TTree* InTree, const TG4Sol_Event& ev,
                      const std::vector<TClonesArray*>& hits)
{
  if(In.status == ReturnRes::EndRun)
    {
      Out.status  = In.status;
      Out.idEvent = In.BeginEvent;
      return;
    }

  auto idEvent = In.BeginEvent;

  InTree->GetEntry(idEvent);
  // ReturnRes::InfoM return_build = (*dynamic_cast<TBuildDetectorLayerPlaneDAF_ZMQ*>(process))(std::ref(ev),
  // std::ref(hits), std::ref(Out));
  ReturnRes::InfoM return_build = (*process)(std::ref(ev), std::ref(hits), std::ref(Out));

  logging->info("f_builder> Out size P{} H{}", Out.DumpParticles.size(), Out.DumpHits.size());
  Out.idEvent = idEvent;
  Out.status  = return_build;
}

struct c_det_build
{
  zmq::context_t& context;
  const std::shared_ptr<spdlog::logger>& logging;
  const std::vector<TDataBuilder*>& list_builder;

  c_det_build(zmq::context_t& c, const std::shared_ptr<spdlog::logger>& l, const std::vector<TDataBuilder*>& b)
      : context(c), logging(l), list_builder(b)
  {
  }

  void operator()(int idTh, TTree* InTree, const TG4Sol_Event& ev, const std::vector<TClonesArray*>& hits,
                  const std::string& addr_IE, const std::string& addr_FB, const std::string& addr_C,
                  const std::string& addr_M)

  {
    zmq::socket_t sock_init_R(context, zmq::socket_type::pull);
    zmq::socket_t sock_init_Rmissing(context, zmq::socket_type::pull);
    zmq::socket_t sockOut(context, zmq::socket_type::push);
    zmq::socket_t control(context, zmq::socket_type::sub);

    static const std::string tag_stats("stats_det_build");
    zmq::socket_t monitor(context, zmq::socket_type::pub);
    
    try
      {
        // std::this_thread::sleep_for(std::chrono::milliseconds(100*idTh));

        sock_init_R.connect(addr_IE.c_str()); //("tcp://127.0.0.1:10001");
                                              // std::this_thread::sleep_for(std::chrono::milliseconds(10*idTh));
        std::string addr_missing(addr_IE);
        addr_missing += "9"; //"Missing";

	sock_init_Rmissing.connect(addr_missing.c_str());

	sockOut.connect(addr_FB.c_str());

	control.connect(addr_C.c_str());
        //control.setsockopt(ZMQ_SUBSCRIBE, "", 0);
	control.set(zmq::sockopt::subscribe, "");
	
	monitor.connect(addr_M.c_str());
	
      }
    catch(zmq::error_t& e)
      {
        logging->error("bind ThQ0 {} : {}", idTh, e.what());
        exit(-1);
      }

    auto f_builder_bind   = std::bind(f_builder, std::placeholders::_1, std::placeholders::_2,
                                    dynamic_cast<TBuildDetectorLayerPlaneDAF_ZMQ*>(list_builder[idTh]), logging, InTree,
                                    std::ref(ev), std::ref(hits));
    using BuilderBindType = decltype(f_builder_bind);

    std::this_thread::sleep_for(std::chrono::milliseconds(200));

    StateMachineThreadIn2Out<EventStatus, DataBuilderOut, BuilderBindType> processB(
        sock_init_R, sock_init_Rmissing, control, sockOut, monitor, tag_stats, f_builder_bind, InOut, "Console",
        "Builder", idTh);
    processB.Run();
  }
};

using ProcessFitter = TDataProcess<ZMQ::DataBuilderOut, ZMQ::DataFitterOut>;

inline void f_RK_daf(ZMQ::DataBuilderOut& In, ZMQ::DataFitterOut& Out, ProcessFitter* process)
{
  if(In.status == ReturnRes::EndRun)
    {
      Out.status  = In.status;
      Out.idEvent = In.idEvent;
      return;
    }
  if(In.status != ReturnRes::Fine)
    {
      Out.status       = In.status;
      Out.idEvent      = In.idEvent;
      Out.previousStep = std::move(In);
      return;
    }

  ReturnRes::InfoM return_process = (*process)(std::ref(In), &Out);
  Out.status                      = return_process;
}

struct c_RK
{

  zmq::context_t& context;
  const std::shared_ptr<spdlog::logger>& logging;

  const std::vector<ProcessFitter*>& list_processMC_RK;

  c_RK(zmq::context_t& c, const std::shared_ptr<spdlog::logger>& l, const std::vector<ProcessFitter*>& r)
      : context(c), logging(l), list_processMC_RK(r)
  {
  }

  void operator()(int idTh, const std::string& addr_BF, const std::string& addr_FF, const std::string addr_C,
                  const std::string& addr_M)
  {
    zmq::socket_t sockQ0_R(context, zmq::socket_type::pull);
    zmq::socket_t sockQout(context, zmq::socket_type::push);
    zmq::socket_t controlQ1(context, zmq::socket_type::sub);

    static const std::string tag_stats("stats_RK");
    zmq::socket_t monitor(context, zmq::socket_type::pub);

    try
      {
        // std::this_thread::sleep_for(std::chrono::milliseconds(100*idTh));

        sockQ0_R.connect(addr_BF.c_str()); //("tcp://127.0.0.1:10001");
        // std::this_thread::sleep_for(std::chrono::milliseconds(10*idTh));
        sockQout.connect(addr_FF.c_str());

        controlQ1.connect(addr_C.c_str());
        //controlQ1.setsockopt(ZMQ_SUBSCRIBE, "", 0);
	controlQ1.set(zmq::sockopt::subscribe, "");

        monitor.connect(addr_M.c_str());
        
      }
    catch(zmq::error_t& e)
      {
        logging->error("bind ThQ1 {} : {}", idTh, e.what());
        exit(-1);
      }

    list_processMC_RK[idTh]->InitMT();

    auto f_RK_bind       = std::bind(f_RK_daf, std::placeholders::_1, std::placeholders::_2, list_processMC_RK[idTh]);
    using FitterBindType = decltype(f_RK_bind);

    std::this_thread::sleep_for(std::chrono::milliseconds(200));

    StateMachineThreadInOut<DataBuilderOut, DataFitterOut, FitterBindType> processRK(
        sockQ0_R, controlQ1, sockQout, monitor, tag_stats, f_RK_bind, InOut, "Console", "RK", idTh);
    processRK.Run();
  }
};

inline void f_merger(TMergerOutput_ZMQ* process, ZMQ::DataFitterOut& In, ZMQ::EventStatus& Out, TTree* OutTree,
                     MCAnaEventG4Sol* OutEvent)
{
  if(In.status == ReturnRes::EndRun)
    {
      Out.BeginEvent = In.idEvent;
      Out.EndEvent   = In.idEvent;
      Out.status     = In.status;
      return;
    }

  Out.BeginEvent = In.idEvent;
  Out.EndEvent   = In.idEvent;
  // ReturnRes::InfoM return_merger = (*dynamic_cast<TMergerOutput_ZMQ*>(process))(std::ref(In), OutEvent);
  ReturnRes::InfoM return_merger = (*process)(std::ref(In), OutEvent);

  OutTree->Fill();
  OutEvent->Clear();

  Out.status = return_merger;
}

struct c_Output
{
  zmq::context_t& context;
  const std::shared_ptr<spdlog::logger>& logging;
  const std::vector<TDataMerger*>& list_merger;

  c_Output(zmq::context_t& c, const std::shared_ptr<spdlog::logger>& l, const std::vector<TDataMerger*>& b)
      : context(c), logging(l), list_merger(b)
  {
  }

  void operator()(int idTh, TTree* OutTree, MCAnaEventG4Sol* OutEvent, const std::string& addr_BM,
                  const std::string& addr_FM, const std::string& addr_C, const std::string& addr_M)
  {
    logging->debug("MT3> start thread Output");

    zmq::socket_t sockQ1_R(context, zmq::socket_type::pull);
    zmq::socket_t sockQout(context, zmq::socket_type::push);
    zmq::socket_t controlQ(context, zmq::socket_type::sub);
    static const std::string tag_stats("stats_merger");
    zmq::socket_t monitor(context, zmq::socket_type::pub);
    try
      {
        // std::this_thread::sleep_for(std::chrono::milliseconds(100*idTh));

        sockQ1_R.connect(addr_BM.c_str()); //("tcp://127.0.0.1:10001");
        // std::this_thread::sleep_for(std::chrono::milliseconds(10*idTh));

        sockQout.connect(addr_FM.c_str());
        controlQ.connect(addr_C.c_str());
        //controlQ.setsockopt(ZMQ_SUBSCRIBE, "", 0);
	controlQ.set(zmq::sockopt::subscribe, "");

        monitor.connect(addr_M.c_str());
        
      }
    catch(zmq::error_t& e)
      {
        logging->error("bind ThQ2 {} : {}", idTh, e.what());
        exit(-1);
      }
    std::this_thread::sleep_for(std::chrono::milliseconds(200));

    auto f_merger_bind = std::bind(f_merger, dynamic_cast<TMergerOutput_ZMQ*>(list_merger[idTh]), std::placeholders::_1,
                                   std::placeholders::_2, OutTree, OutEvent);
    using MergerBindType = decltype(f_merger_bind);

    StateMachineThreadInOut<DataFitterOut, EventStatus, MergerBindType> processM(
        sockQ1_R, controlQ, sockQout, monitor, tag_stats, f_merger_bind, InOut, "Console", "Merger",
        idTh); //, sizeQ2, sizeQend);
    processM.Run();
  }
};

#endif
