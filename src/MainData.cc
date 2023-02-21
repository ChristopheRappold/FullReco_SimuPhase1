#include "Ana_Event/Ana_WasaEvent.hh"
#include "Ana_Hist.hh"
#include "Debug.hh"
#include "FullRecoConfig.hh"
#include "FullRecoTask.hh"
//#include "FullRecoTaskMT.hh"
//#include "FullRecoTaskZMQ.hh"
#include "KFParticle.h"


//#include "TApplication.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TEntryList.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TKey.h"
//#include "TMacro.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TTree.h"
#include "TTreePerfStats.h"
#include "TROOT.h"

#include <iostream>
#include <string>
#include <vector>

#ifdef ROOT6
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#endif

#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

using namespace std;

int main(int argc, char** argv)
{
  std::cout << " Config : ";
  FullRecoConfig config(argc, argv);
  if(config.ProperConf() != 0)
    return -1;
  std::cout << config.CheckConfig();

  try
    {
      int Nb_CPU           = config.Get<int>("Nb_CPU");
      int Nb_event         = config.Get<int>("Nb_Event_Cfg");
      int Nb_fraction      = config.Get<int>("Nb_Fraction");
      int Start            = config.Get<int>("Start");
      int Log_lvl          = config.Get<int>("Log_Lvl");
      std::string nameGeo  = config.Get<std::string>("Geo");
      std::string name_in  = config.Get<std::string>("Input_Namefile");
      std::string name_out = config.Get<std::string>("Output_Namefile");
      bool Restarting = config.IsAvailable("Task_ReStart") ? config.Get<bool>("Task_ReStart") : false;

      if(config.IsAvailable("LocalOutput_Suffix"))
        {
          std::string file_base_name = name_out.substr(0,name_out.find_last_of('.'));
          file_base_name += config.Get<std::string>("LocalOutput_Suffix");
          file_base_name += ".root";
          name_out = file_base_name;
        }

      bool MT = config.IsAvailable("MultiThreading");
      bool ZMQ = config.IsAvailable("ZeroMQ");
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

      KFParticle test_particle;
      if(test_particle.IsHomogeneous() == 1)
        ConsoleLogger->info("KFParticle Magnetic field: Homogeneous");
      else if(test_particle.IsHomogeneous() == 10)
        ConsoleLogger->info("KFParticle Magnetic field: NonHomogeneous");

      std::string file_base_name = name_out.substr(name_out.find_last_of('/') + 1);
      ConsoleLogger->info("name in : {} / name out : {}", name_in, name_out);

      /*******************************************************************/

      TStopwatch timer;
      timer.Start(kFALSE);

      // HypHiMC_output_TREE ***********************
      Long64_t total_nentries = 0;
      // Long64_t total_nentries_over_files=0;
      Long64_t Start_event = 0;
      Long64_t Stop_event = 0;

      TGeoManager::Import(nameGeo.c_str()); // For Kalman_DAF

      TDatabasePDG::Instance()->AddParticle("deuteron", "deuteron", 1.875477 /*1.875613*/, kTRUE, 0., 1. * 3., "Ions", 10000);
      TDatabasePDG::Instance()->AddParticle("triton", "triton", 2.807971/*2.80892*/ /*2.80925*/, kTRUE, 0., 1. * 3., "Ions", 10001);
      TDatabasePDG::Instance()->AddParticle("alpha", "alpha", 3.726954 /*3.72738*/ /*3.727417*/, kTRUE, 0., 2. * 3., "Ions", 10002);
      TDatabasePDG::Instance()->AddParticle("He3", "He3", 2.80746 /*2.80839*/ /*2.80923*/, kTRUE, 0., 2. * 3., "Ions", 10003);
      TDatabasePDG::Instance()->AddParticle("Li6", "Li6", 5.601432 /*5.60152*/ /*5.6015194*/, kTRUE, 0., 3. * 3., "Ions", 10004);
      TDatabasePDG::Instance()->AddParticle("H3L", "H3L", 2.99114, kFALSE, 0., 1. * 3., "Ions", 20001);
      TDatabasePDG::Instance()->AddParticle("H4L", "H4L", 3.9225, kFALSE, 0., 1. * 3., "Ions", 20002);
      TDatabasePDG::Instance()->AddParticle("He5L", "He5L", 4.8399, kFALSE, 0., 2. * 3., "Ions", 20003);

      TFile* offile = new TFile(name_out.c_str(), "RECREATE");

      TTree* Anatree          = 0;
      Ana_WasaEvent* Anaevent = 0;
      Anatree                 = new TTree("T", "WASAUnpackTree");
      // Anatree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
      Anatree->SetAutoSave(100000000); // autosave when 100 Mbyte written
      // Anatree->SetAutoSave(50000000); // autosave when 50 Mbyte written
      // Anatree->SetMaxTreeSize(50000000); //1 GB
      Anatree->SetAutoFlush();
      Anaevent = new Ana_WasaEvent();
      Ana_WasaEvent::Class()->IgnoreTObjectStreamer();
      TBranch* branch = Anatree->Branch("Ana_WasaEvent", &Anaevent, 32000, 2);
      branch->SetAutoDelete(kFALSE);
      Anatree->BranchRef();
      Anatree->SetDirectory(offile);

      long long int nb = 0;

      TFile* input_file = new TFile(name_in.c_str());

      TTree* InTree = nullptr;
      if(Restarting == false)
        InTree        = input_file->Get<TTree>("WASAUnpackTree");
      else
        InTree        = input_file->Get<TTree>("T");

      assert(InTree != nullptr);

      DataSimExp InputPar{nullptr, nullptr, nullptr};
      //InputPar.simexpMetadata  = (AnaEvent_Metadata*)(input_file->Get("EventMetadata"));
/*
      if(InputPar.simexpMetadata == nullptr)
        {
          ConsoleLogger->error("E> Load AnaEvent_Metadata not possible ! {}", fmt::ptr(input_file->Get("EventMetadata")));
          return -1;
        }
*/
      EventWASAUnpack fEvent ;
      Ana_WasaEvent* fEventRe = new Ana_WasaEvent();

      if(Restarting == false)
        {
          InTree->SetBranchAddress("lmd_header",&fEvent.header);
          InTree->SetBranchAddress("s4tq"      ,&fEvent.s4tq);
          InTree->SetBranchAddress("s4mwdc"    ,&fEvent.s4mwdc);
          InTree->SetBranchAddress("s4wfd"     ,&fEvent.s4wfd);
          InTree->SetBranchAddress("s2tq1"     ,&fEvent.s2tq1);
          InTree->SetBranchAddress("s2mdc"     ,&fEvent.s2mdc);
          InTree->SetBranchAddress("s2wfd123"  ,&fEvent.s2wfd123);
          InTree->SetBranchAddress("s2tq2"     ,&fEvent.s2tq2);
          InTree->SetBranchAddress("s2wfd45"   ,&fEvent.s2wfd45);
          InTree->SetBranchAddress("s2fiber"   ,&fEvent.s2fiber);
          InTree->SetBranchAddress("s2csi"     ,&fEvent.s2csi);
          InTree->SetBranchAddress("frstpc"    ,&fEvent.frstpc);
        }
      else
      	InTree->SetBranchAddress("Ana_WasaEvent", &fEventRe);


      Long64_t total_nentries_over_files = InTree->GetEntries();

      Start_event = Start;
      Stop_event  = Nb_event != -1 ? Nb_event + Start : total_nentries_over_files;

      if(Stop_event > total_nentries_over_files)
        Stop_event = total_nentries_over_files;
      if(Start_event > total_nentries_over_files)
        Start_event = total_nentries_over_files;


      config.Add("Nb_Event", Stop_event - Start_event);
      config.Add("Start_Event", Start_event);
      config.Add("Stop_Event", Stop_event);

      Ana_Hist ListHisto(true /*Daf*/, false /*Oldvertex*/, false /*DCproject*/, true /*Finding*/, true /*Riemann*/, false /*Hough*/,
                  false /*Simu*/, true /*Builder*/, true /*PrimVtx*/, false /*PrimVtx_Si*/, true /*DecayVtx*/, true /*FragmentFinder*/,
                  true /*WASAFinder*/);

      AnaEvent_Metadata metadata;

      if(MT == false && ZMQ == false)
        {
          FullRecoTask<Ana_WasaEvent> ReconstructionTask(config, InputPar);
          
          ReconstructionTask.AttachHisto(&ListHisto);
          ConsoleLogger->info("Init done");
          //-------------------------------------------------------
          // EVENT LOOP nentries
          TTimeStamp stamp;
          
          ConsoleLogger->info("Start = {} Stop= {}", Start_event, Stop_event);
          ConsoleLogger->info("---------------------------------------------");
          
          Long64_t for_display = 0;
          Long64_t timing      = 0;
          // reader.SetEntriesRange(Start_event,Stop_event);
          
          for(Long64_t i = Start_event; i < Stop_event; ++i)
            {
              InTree->GetEntry(i);
              Long64_t iEvent = i; // reader.GetCurrentEntry();

              ++total_nentries;
              Anaevent->Clear();

              if(iEvent % 10000 == 0)
                ConsoleLogger->info("Processing Event# {} / {} | {} %", iEvent, Stop_event, (double)iEvent / (double)(Stop_event - Start_event) * 100);

              if((int)((double)iEvent / (double)(Stop_event - Start_event) * 4) == timing)
                {
                  ConsoleLogger->info("Progress : {} %",
                  (int)((double)iEvent / (double)(Stop_event - Start_event) * 100.));
                  ++timing;
                }

              int toStop = 0;
              if(Restarting == false)
                toStop = ReconstructionTask.EventLoop(fEvent, Anaevent);
              else
                toStop = ReconstructionTask.EventLoop(fEventRe, Anaevent);

              nb += Anatree->Fill();
              ++for_display;
              Anaevent->Clear();
            }

          ReconstructionTask.SetEventMetadata(metadata);
        }

      else
        ConsoleLogger->warn("W> no run : MT? {} ZMQ? {}",MT,ZMQ);

      offile = Anatree->GetCurrentFile();
      offile->cd();
      Anatree->Write();
      ConsoleLogger->info("Tree Filled");
      ListHisto.Write(offile);

      ConsoleLogger->info("Histo Written");

      TObjArray* additionalPDG = new TObjArray;
      additionalPDG->SetName("additionalPDG");
      int id_pdg  = 10004 + 1;
      bool finish = false;
      while(finish == false)
        {
          TParticlePDG* PDG_particle = TDatabasePDG::Instance()->GetParticle(id_pdg);
          if(PDG_particle != nullptr)
            {
              additionalPDG->Add(PDG_particle);
              ++id_pdg;
            }
          else
            finish = true;
        }
      id_pdg = 20003 + 1;
      while(finish == false)
        {
          TParticlePDG* PDG_particle = TDatabasePDG::Instance()->GetParticle(id_pdg);
          if(PDG_particle != nullptr)
            {
              additionalPDG->Add(PDG_particle);
              ++id_pdg;
            }
          else
            finish = true;
        }
      offile->WriteObjectAny(additionalPDG, additionalPDG->Class(), "additionalPDG");
      ConsoleLogger->info("Additional PDG data written");

      offile->WriteTObject(&metadata,"EventMetadata");
      ConsoleLogger->info("Event metadata written");

      offile->Close();
      ConsoleLogger->info("Histo File Closed");

#ifdef TREEPERF
      ps->SaveAs("perfstat_2.root");
#endif

      delete offile;
      offile = NULL;
      ConsoleLogger->info("Offile deleted");
      // delete offile_hist;offile_hist=NULL;
      // std::cout<<"Offile Histo deleted"<<std::endl;

      delete Anaevent;
      Anaevent = NULL;
      ConsoleLogger->info("Anaevent deleted");

      
      additionalPDG->Delete();
      additionalPDG = nullptr;

#ifdef TREEPERF
      delete ps;
      ps = 0;
#endif
      /***************************************************************/
      timer.Stop();

      Float_t mbytes = 0.000001 * nb;
      Double_t rtime = timer.RealTime();
      Double_t ctime = timer.CpuTime();
      ConsoleLogger->info("");
      ConsoleLogger->info("{} events and {} bytes processed.", total_nentries, nb);
      ConsoleLogger->info("RealTime= {} seconds, CpuTime= {} seconds", rtime, ctime);
      ConsoleLogger->info("You read {} Mbytes/Realtime seconds", mbytes / rtime);
      ConsoleLogger->info("You read {} Mbytes/Cputime seconds", mbytes / ctime);

      return 0;
    }
  catch(const spdlog::spdlog_ex& ex)
    {
      std::cout << "Log initialization failed: " << ex.what() << std::endl;
      std::exit(-1);
    }
}
