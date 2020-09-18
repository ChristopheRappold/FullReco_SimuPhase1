//#include "TApplication.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TROOT.h"
#include "TH1I.h"
#include "TH1F.h"
#include <TFile.h>
#include <TMemFile.h>
#include <THttpServer.h>

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include "FullRecoConfig.hh"

#include "zmq.hpp"
#include "zmq_addon.hpp"
#include "msgpack.hpp"

using namespace std;

template <typename T>
void subscribe(zmq::socket_t& socket, std::string& tag, T& data)
{

  zmq::multipart_t all(socket);

  tag = std::string(static_cast<const char*>(all.at(0).data()), all.at(0).size());

  msgpack::unpacked unpacked_body;
  msgpack::object_handle unpack_body = msgpack::unpack(static_cast<const char*>(all.at(1).data()), all.at(1).size());
  msgpack::object obj = unpack_body.get();
  obj.convert(data);
}


int main(int argc, char** argv)
{
  std::cout << " Config : ";
  FullRecoConfig config(argc, argv);
  if(config.ProperConf() != 0)
    return -1;
  std::cout << config.CheckConfig();

  try
    {
      bool MT = config.IsAvailable("MultiThreading");
      bool ZMQ = config.IsAvailable("ZeroMQ");
      int Log_lvl          = config.Get<int>("Log_Lvl");
      std::string addr_monitor = config.IsAvailable("Addr_Monitor") ? config.Get<std::string>("Addr_Monitor") : "tcp://127.0.0.1:9876";
      std::string name_in  = config.Get<std::string>("Input_Namefile");
      
      if(MT || ZMQ)
	{
	  //ROOT::EnableThreadSafety();
	  //ROOT::EnableImplicitMT();
	}

      auto sink1 = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
      sink1->set_pattern("[%^%l%$] %v");
      auto ConsoleLogger = std::make_shared<spdlog::logger>("Console", sink1);
      spdlog::register_logger(ConsoleLogger);
      // auto console = spdlog::stdout_color_mt("Console");
      switch(Log_lvl)
        {
        case -1:
          ConsoleLogger->set_level(spdlog::level::off);
          break;
        case 0:
          ConsoleLogger->set_level(spdlog::level::warn);
          break;
        case 1:
          ConsoleLogger->set_level(spdlog::level::info);
          break;
        case 2:
          ConsoleLogger->set_level(spdlog::level::debug);
          break;
        default:
          ConsoleLogger->set_level(spdlog::level::info);
          break;
        }

      ConsoleLogger->info("Multi-threading :{}", MT);

      THttpServer* serv = new THttpServer(Form("http:8080?top=%s", "stats"));
      serv->SetReadOnly(kFALSE);

      ConsoleLogger->info("Add Monitor :{}",addr_monitor);
      //TFile* file = new TMemFile(name_in.c_str(),"RECREATE","Stats of the run"); 
      TH1I* h1 = new TH1I("MonitorStats","MonitorStats",10,0,10);
      h1->SetDirectory(0);
      serv->Register("/",h1);

      struct Hist {
	int previous_bin;
	double previous_content;
	TH1F* h_duration;
	TH1F* h_scaler;
	
      };
      
      std::unordered_map<std::string, Hist> hists;
      auto Fill_Hists = [&hists,&serv,&ConsoleLogger](const std::string& idH, auto bin, auto content) {
	
	auto it_h = hists.find(idH);
	if(it_h == hists.end())
	  {
	    ConsoleLogger->info("new Histo : {} / {} {}",idH, bin, content);

	    std::string nameD(idH);
	    static const std::string duration("_duration");
	    nameD += duration;

	    TH1F* h = new TH1F(nameD.c_str(), nameD.c_str(), 1000,0,1000);

	    h->SetDirectory(0);
	    h->SetOption("hist");
	    h->SetCanExtend(TH1::kXaxis); 

	    std::string name_dir("/");
	    name_dir += idH;
	    serv->Register(name_dir.c_str(),h);

	    h->Fill(bin, content);

	    
	    std::string nameS(idH);
	    static const std::string scaler("_scaler");
	    nameS += scaler;
	    TH1F* h2 = new TH1F(nameS.c_str(), nameS.c_str(), 1000,0,1000);
	    
	    h2->SetDirectory(0);
	    h2->SetOption("hist");
	    h2->SetCanExtend(TH1::kXaxis); 

	    serv->Register(name_dir.c_str(),h2);

	    auto freq = bin/static_cast<double>(content);
	    int time_init = 0;
	    int time_final = content/10;

	    for(int i=time_init;i<=time_final;++i)
	      h2->SetBinContent(i,freq/10);
	    
	    Hist hist = {bin, content, h, h2};
	    
	    auto histo = hists.insert(std::make_pair(idH, hist));
	  }
	else
	  {
	    auto& hist = it_h->second;
	    auto new_bin = bin;
	    auto new_content = content;
	    
	    hist.h_duration->Fill(new_bin,new_content);

	    auto freq = (bin-hist.previous_bin)/static_cast<double>(content);
	    int time_init = hist.previous_content/10;
	    int time_final = (hist.previous_content+content)/10;
	    
	    for(int i=time_init;i<time_final;++i)
	      hist.h_scaler->SetBinContent(i,freq/10); 
	    
	    hist.previous_bin = new_bin;
	    hist.previous_content += content;	    
	  }
	//file->Write();
      };

      auto f_zmq = [&]() {

      
	zmq::context_t context;
	zmq::socket_t socket(context,  zmq::socket_type::sub);
	socket.bind(addr_monitor.c_str());
	const std::array<std::string,5> tags = {"stats_nEvent","stats_fEvent","stats_det_build","stats_RK","stats_merger"};
	for(auto& tag : tags) 
	  socket.set(zmq::sockopt::subscribe, tag);//.c_str(), tag.size());

	//socket.setsockopt(ZMQ_LINGER, 1000);
	socket.set(zmq::sockopt::linger, 1000);

	
      std::string dash_char("_");
      //int count = 0;
      while(1)
	{
	  std::string tag_data;
	  std::tuple<int, int, int, double> data;
	  subscribe(socket,tag_data,data);

	  ConsoleLogger->info("stats> {} : {}/{} | {} in {} second",tag_data, std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data));

	  h1->Fill(tag_data.c_str(),1.);
	  //file->Write();
	  tag_data += std::to_string(std::get<0>(data));
	  tag_data += dash_char;
	  tag_data += std::to_string(std::get<1>(data));
	  Fill_Hists(tag_data,std::get<2>(data),std::get<3>(data));

	  //++count;
	  // if(count==20)
	  //   {
	  //     //file->Write();
	  //     count = 0;
	  //   }

	}

      };

      //serv->SetJSROOT("http://jsroot.gsi.de/latest/");

      // delete AnaOutTree;
      std::thread t1(f_zmq);
      t1.detach();
 
      const Long_t kUPDATE = 300;
      Long_t cnt = 0;
      
      while(1)
	{
	  if (cnt++ % kUPDATE == 0)
	    if (gSystem->ProcessEvents())
	      break;
	}
      
      return 0;
    }
  catch(const spdlog::spdlog_ex& ex)
    {
      std::cout << "Log initialization failed: " << ex.what() << std::endl;
      std::exit(-1);
    }
}
