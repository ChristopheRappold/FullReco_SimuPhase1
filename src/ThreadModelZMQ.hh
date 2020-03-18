#ifndef TTHREADMODEL_MT
#define TTHREADMODEL_MT

#include "msgpack.hpp"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"
#include "zmq.hpp"
#include "zmq_addon.hpp"

#include <algorithm>
#include <chrono>
#include <type_traits>
#include <iostream>
#include <thread>
#include <sys/types.h>
#include <unistd.h>

#include "ReturnRes.hh"

enum States : int
{
  init = 0,
  process,
  stats,
  finishing,
  last,
  end,
};

enum typeIO : int
{
  OutOnly = 0,
  InOut,
  InOnly,
};

template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())>
{
};
// For generic types, directly use the result of the signature of its 'operator()'

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType (ClassType::*)(Args...) const>
// we specialize for pointers to member function
{
  enum
  {
    arity = sizeof...(Args)
  };
  // arity is the number of arguments.

  typedef ReturnType result_type;

  template <size_t i>
  struct arg
  {
    typedef typename std::tuple_element<i, std::tuple<Args...> >::type type;
    // the i-th argument is equivalent to the i-th tuple element of a tuple
    // composed of those arguments.
  };
};

template <typename F>
using In_arg = std::decay_t<typename function_traits<F>::template arg<0>::type>;
template <typename F>
using Out_arg = std::decay_t<typename function_traits<F>::template arg<1>::type>;


template <typename T, typename = int>
struct HasEndEvent : std::false_type { };
template <typename T>
struct HasEndEvent <T, decltype((void) T::EndEvent, 0)> : std::true_type { };

template <typename T, typename = int>
struct HasIdEvent : std::false_type { };
template <typename T>
struct HasIdEvent <T, decltype((void) T::idEvent, 0)> : std::true_type { };

template<typename In>
std::tuple<Long64_t,Long64_t> RangeEvent(const In& Status, typename std::enable_if<HasEndEvent<In>::value>::type * = nullptr )
{
  //std::cout<<"RangeEvent B/E:"<<Status.BeginEvent<<" "<<Status.EndEvent<<"\n";
  return std::make_tuple(Status.BeginEvent,Status.EndEvent);
}
template<typename In>
std::tuple<Long64_t,Long64_t> RangeEvent(const In& Status, typename std::enable_if<HasIdEvent<In>::value>::type * = nullptr )
{
  //std::cout<<"RangeEvent id:"<<Status.idEvent<<"\n";
  return std::make_tuple(Status.idEvent,Status.idEvent+1);
}

template<typename In>
void SetEvent(Long64_t i, In& Status, typename std::enable_if<HasEndEvent<In>::value>::type * = nullptr )
{
  Status.BeginEvent = i;
  //Status.EndEvent = i;
}
template<typename In>
void SetEvent(Long64_t i, In& Status, typename std::enable_if<HasIdEvent<In>::value>::type * = nullptr )
{
  Status.idEvent = i;
}

template <typename T>
void publish(zmq::socket_t& socket, const std::string& tag, const T& data)
{
  msgpack::sbuffer packed;
  msgpack::pack(&packed, data);

  zmq::multipart_t all_data;
  all_data.addstr(tag);
  all_data.addmem(packed.data(),packed.size());

  all_data.send(socket);
}



template <typename In, typename Out, typename Functor, int Nstats = 100>
// template<typename Functor>
class StateMachineThread
{
public:
  // using In = In_arg<Functor>;
  // using Out = Out_arg<Functor>;
  zmq::socket_t& monitor;
  const std::string& tag_stats;

  Functor& processor;
  typeIO type;
  std::shared_ptr<spdlog::logger> logging;
  std::string name;
  int idTh;

  States state = init;
  int count    = 0;
  int previousStats = 0;
  
public:
  virtual void Send(const Out& out, const std::string& addlog){};
  virtual void Receive(In& in, const std::string& addlog){};
  StateMachineThread() = delete;
  StateMachineThread( zmq::socket_t& _sockMonitor, const std::string& tagM,//, zmq::socket_t& _sockOut,
      Functor& _f, typeIO _type, const std::string& nameLog, const std::string& _name, int _idTh)
    : monitor(_sockMonitor), tag_stats(tagM), processor(_f), type(_type), logging(spdlog::get(nameLog.c_str())), name(_name), idTh(_idTh)
  {
  }
  
  void Run()
  {
    auto pidSys = ::getpid();//std::hash<std::thread::id>()( std::this_thread::get_id() );

    while(1)
      {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
	auto time_ref_start = std::chrono::high_resolution_clock::now();
	
        switch(state)
          {
          case init:
            {
              logging->info("MT{}> start thread {} / id {}", name, name, idTh);
              state = process;
            }
          case process:
            {
              In Status;
              if(type == InOnly || type == InOut)
                Receive(Status, "receiving");

              if(Status.status == ReturnRes::CloseRun)
                {
                  logging->info("MT{}> got status CloseRun ! | id:{} ", name, idTh);
                  state = end;
                  break;
                }
	      auto [first, last] = RangeEvent(Status);
	      
	      logging->info("MT{}> event#{} / [B{},E{},S{}] / {}", name, first, first, last, Status.status ,
                            idTh); //,isR->size,isR->untruncated_size);

	      for(auto iter = first; iter<last; ++iter)
		{
		  SetEvent(iter,Status);
		  logging->debug("MT{}> In loop #{} / [B{},E{},S{}] / {}", name, iter, first, last, Status.status ,idTh);
		  Out StatusO;
		  processor(Status, StatusO);
		  
		  if(type == OutOnly || type == InOut)
		    Send(StatusO, "send");

		  ++count;
		}
              if(count == Nstats)
                {
                  state = stats;
                  count = 0;
                }
              break;
            }
          case stats:
            {
	      auto time_stat = std::chrono::high_resolution_clock::now();
	      std::chrono::duration<double> time_diff = time_stat - time_ref_start; 
	      
	      std::tuple<int,int,int,double> time_data = std::make_tuple(pidSys, idTh, Nstats*(previousStats+1), std::chrono::duration_cast<std::chrono::nanoseconds>(time_diff).count());
	      publish(monitor,tag_stats,time_data);
	      ++previousStats;
	      time_ref_start = time_stat;
	      state = process;
              break;
            }
          case finishing:
            {
              logging->info("MT{}> finishing ! | id:{}", name, idTh);
              In Status;
              if(type == InOnly || type == InOut)
                Receive(Status, "receiving");

              if(Status.status == ReturnRes::CloseRun)
                {
                  logging->info("MT{}> got status CloseRun ! | id:{} tempS:{}", name, idTh); //, tempS);
                  state = end;
                  break;
                }

              logging->info("MT{}> got message ! id:{} {}", name, idTh, Status.status);
              break;
            }
          case end:
            logging->info("MT{}> End ! Th#{}", name, idTh);
            return;

          case last:
            {
              logging->info("MT{}> finishing ! | push Q next finished status Th#{}", name, idTh); //, tempS);
              state = end;
            }
          }
      }
  }
};

template <typename In, typename Out, typename Functor, int Nstats = 100>
// template<typename Functor>
class StateMachineThreadInOut : public StateMachineThread<In, Out, Functor, Nstats>
{
public:
  // using In = In_arg<Functor>;
  // using Out = Out_arg<Functor>;

  zmq::socket_t& sockIn;
  zmq::socket_t& control;
  zmq::socket_t& sockOut;
  zmq::active_poller_t pollSockets;

  bool message_received1 = false;
  bool message_received2 = false;
 
public:
  StateMachineThreadInOut() = delete;
  StateMachineThreadInOut(zmq::socket_t& _sockIn, zmq::socket_t& _control, zmq::socket_t& _sockOut, zmq::socket_t& _sockMonitor, const std::string& tag,
			  Functor& _f, typeIO _type, const std::string& nameLog, const std::string& _name, int _idTh)
    : StateMachineThread<In, Out, Functor>(_sockMonitor, tag, _f, _type, nameLog, _name, _idTh), sockIn(_sockIn), control(_control),
        sockOut(_sockOut), pollSockets()
  {
    bool& m_rcv1 = message_received1;

    zmq::active_poller_t::handler_type handler1 = [&m_rcv1](zmq::event_flags events) {
      if(zmq::event_flags::none != (events & zmq::event_flags::pollin))
        m_rcv1 = true;
    };

    bool& m_rcv2 = message_received2;

    zmq::active_poller_t::handler_type handler2 = [&m_rcv2](zmq::event_flags events) {
      if(zmq::event_flags::none != (events & zmq::event_flags::pollin))
        m_rcv2 = true;
    };

    pollSockets.add(sockIn, zmq::event_flags::pollin, handler1);
    pollSockets.add(control, zmq::event_flags::pollin, handler2);
  }

  void Receive(In& Status, const std::string& addlog)
  {
    try
      {
        int rs = pollSockets.wait(std::chrono::milliseconds{-1});
        if(rs > 0)
          {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            if(message_received1)
              {
                zmq::message_t body_msg;
                sockIn.recv(body_msg); //,zmq::recv_flags::none);//bufIn,zmq::recv_flags::dontwait);

                msgpack::object_handle unpacked_body =
                    msgpack::unpack(static_cast<const char*>(body_msg.data()), body_msg.size());
                msgpack::object obj = unpacked_body.get();
                obj.convert(Status);
		this->logging->info("MT{}#{} received message !", this->name, this->idTh);

                message_received1 = false;
              }
            if(message_received2)
              {
                this->logging->info("MT{}#{} Terminate message !", this->name, this->idTh);
                zmq::message_t body_msg;
                control.recv(body_msg);
                const std::string str_msg(static_cast<const char*>(body_msg.data()), body_msg.size());
                static const std::string terminate = "TERMINATE";
                if(str_msg == terminate)
                  {
                    SetEvent(0,Status);
                    Status.status = ReturnRes::CloseRun;
		    
                    this->logging->info("MT{}#{} status break !", this->name, this->idTh);
                  }
                else
                  this->logging->warn("MT{}#{} did not match the msg: {} ", this->name, this->idTh, str_msg);
                message_received2 = false;
              }
          }
      }
    catch(zmq::error_t& e)
      {
        this->logging->error("Error catched in poller MT{}#{} {} ! {} {}", this->name, this->idTh, addlog, e.what(),
                             e.num());
      }
  }
  void Send(const Out& StatusO, const std::string& addlog)
  {
    try
      {
        msgpack::sbuffer sbuf;
        msgpack::pack(sbuf, StatusO);

        zmq::const_buffer body_msg(sbuf.data(), sbuf.size());

        sockOut.send(body_msg, zmq::send_flags::dontwait);
      }
    catch(zmq::error_t& e)
      {
        this->logging->error("Error catched in MT{}#{} {} ! {} {}", this->name, this->idTh, addlog, e.what(), e.num());
      }
  }
};


template <typename In, typename Out, typename Functor, int Nstats = 100>
// template<typename Functor>
class StateMachineThreadIn2Out : public StateMachineThreadInOut<In, Out, Functor, Nstats>
{
  zmq::socket_t& sockIn2;
  bool message_received3 = false;

public:
  StateMachineThreadIn2Out() = delete;
  StateMachineThreadIn2Out(zmq::socket_t& _sockIn, zmq::socket_t& _sockIn2, zmq::socket_t& _control, zmq::socket_t& _sockOut,
			   zmq::socket_t& _sockM, const std::string& tag,
			   Functor& _f, typeIO _type, const std::string& nameLog, const std::string& _name, int _idTh)
    : StateMachineThreadInOut<In, Out, Functor>(_sockIn, _control, _sockOut, _sockM, tag, _f, _type, nameLog, _name, _idTh), sockIn2(_sockIn2)
  {
    // bool& m_rcv1 = this->message_received1;

    // zmq::active_poller_t::handler_type handler1 = [&m_rcv1](zmq::event_flags events) {
    //   if(zmq::event_flags::none != (events & zmq::event_flags::pollin))
    //     m_rcv1 = true;
    // };

    // bool& m_rcv2 = this->message_received2;

    // zmq::active_poller_t::handler_type handler2 = [&m_rcv2](zmq::event_flags events) {
    //   if(zmq::event_flags::none != (events & zmq::event_flags::pollin))
    //     m_rcv2 = true;
    // };
    
    bool& m_rcv3 = message_received3;

    zmq::active_poller_t::handler_type handler3 = [&m_rcv3](zmq::event_flags events) {
      if(zmq::event_flags::none != (events & zmq::event_flags::pollin))
        m_rcv3 = true;
    };

    //this->pollSockets.add(this->sockIn, zmq::event_flags::pollin, handler1);
    //this->pollSockets.add(this->control, zmq::event_flags::pollin, handler2);
    this->pollSockets.add(sockIn2, zmq::event_flags::pollin, handler3);
  }

  void Receive(In& Status, const std::string& addlog)
  {
    try
      {
        int rs = this->pollSockets.wait(std::chrono::milliseconds{-1});
        if(rs > 0)
          {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            if(this->message_received1)
              {
                zmq::message_t body_msg;
                this->sockIn.recv(body_msg); //,zmq::recv_flags::none);//bufIn,zmq::recv_flags::dontwait);

                msgpack::object_handle unpacked_body =
                    msgpack::unpack(static_cast<const char*>(body_msg.data()), body_msg.size());
                msgpack::object obj = unpacked_body.get();
                obj.convert(Status);
		this->logging->info("MT{}#{} received message !", this->name, this->idTh);

                this->message_received1 = false;
              }
            if(this->message_received2)
              {
                this->logging->info("MT{}#{} Terminate message !", this->name, this->idTh);
                zmq::message_t body_msg;
                this->control.recv(body_msg);
                const std::string str_msg(static_cast<const char*>(body_msg.data()), body_msg.size());
                static const std::string terminate = "TERMINATE";
                if(str_msg == terminate)
                  {
                    SetEvent(0,Status);
                    Status.status = ReturnRes::CloseRun;
		    
                    this->logging->info("MT{}#{} status break !", this->name, this->idTh);
                  }
                else
                  this->logging->warn("MT{}#{} did not match the msg: {} ", this->name, this->idTh, str_msg);
                this->message_received2 = false;
              }
	    if(message_received3)
              {
                zmq::message_t body_msg;
                sockIn2.recv(body_msg); //,zmq::recv_flags::none);//bufIn,zmq::recv_flags::dontwait);

                msgpack::object_handle unpacked_body =
		  msgpack::unpack(static_cast<const char*>(body_msg.data()), body_msg.size());
                msgpack::object obj = unpacked_body.get();
                obj.convert(Status);
		this->logging->info("MT{}#{} received message bis !", this->name, this->idTh);

                message_received3 = false;
              }
          }
      }
    catch(zmq::error_t& e)
      {
        this->logging->error("Error catched in poller MT{}#{} {} ! {} {}", this->name, this->idTh, addlog, e.what(),
                             e.num());
      }
  }

  
};


class proxyPullPushNN {

  zmq::context_t& context;
  std::shared_ptr<spdlog::logger> logging;
  
  zmq::socket_t frontend;
  zmq::socket_t backend;
  zmq::socket_t control;

  std::string infoMsg;
public:
  proxyPullPushNN() = delete;
  proxyPullPushNN(zmq::context_t& c, const std::string& addr_frontend, const std::string& addr_backend, const std::string& addr_control,  const std::string& log, const std::string& infoM):
    context(c),logging(spdlog::get(log.c_str())),
    frontend(context, zmq::socket_type::pull),backend(context, zmq::socket_type::push),control(context, zmq::socket_type::sub),infoMsg(infoM)

  {
    frontend.bind(addr_frontend.c_str());
    backend.bind(addr_backend.c_str());
    control.connect(addr_control.c_str());
    control.setsockopt(ZMQ_SUBSCRIBE, "", 0);
  }

  void Run()
  {
    logging->info("Mproxy> proxy {} setting :",infoMsg);
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    try
      {
	zmq::proxy_steerable(frontend,backend,zmq::socket_ref(),control);
      }
    catch(zmq::error_t& e)
      {
	logging->error("Mproxy> {} : {}",infoMsg,e.what());
      }
    //logging->info("Mproxy> -> done !");
  }
};

#endif
