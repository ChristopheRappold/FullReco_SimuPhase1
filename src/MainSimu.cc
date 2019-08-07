#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cstdlib>

#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TChain.h"
#include "TEntryList.h"


#include "TStopwatch.h"

#include "Ana_Hist.hh"
#include "FullRecoTask.hh"

#include "TTimeStamp.h"

#include "TGeoManager.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TDatabasePDG.h"

#include "Debug.hh"

#ifdef ROOT6
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#endif

#include "TTreePerfStats.h"

#include "TMacro.h"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/ostream_sink.h"

using namespace std;
//class spdlog::logger;

static struct option optlong[] = 
  {
    {"help",0,NULL,'h'},
    {"cpu",1,NULL,'c'},
    {"num",1,NULL,'n'},
    {"event",1,NULL,'e'},
    {"start",1,NULL,'s'},
    {"geo",1,NULL,'g'},
    {"log",1,NULL,'l'}
  };

int main(int argc,char** argv)
{
  //int NofLevels = 0;

  if (argc < 3)
    {
      std::cout << "!> Wrong number of parameters!\n";
      std::cout << "!> Example of use:\n";
      std::cout << "!> " << argv[0];
      std::cout << "[-c nb_cpu] [--cpu nb_cpu] [-n fraction] [--num fraction] [-s start_ev] [--start start_ev] [-e nb_event] [--event nb_event] [-l lvllog] [--log lvllog] [-h] OutputFile RootInputFile_withParObj \n";
      std::exit(1);
    }

  int option_char;
  int Nb_CPU=1;
  int Nb_event=-1;
  int Nb_fraction=1;
  int Start=0;
  std::string nameEventList("no");
  bool InputSimuTree = false;
  bool EventGenConf = false;
  TString file_conf;
  bool Mocadi = false;
  bool HypHI = false;
  int Log_lvl = 1;
  std::string nameGeo("./geo/GeoSolenoid.root");
  //int type_hyp = 0;
  while ((option_char = getopt_long (argc,argv,"+hc:n:e:s:g:l:",optlong,NULL)) != EOF)
    switch (option_char)
      {  
      case 'h': std::cerr <<"usage: "<<argv[0]<<" [-c nb_cpu] [--cpu nb_cpu] [-n fraction] [--num fraction] [-s start_ev] [--start start_ev] [-e nb_event] [--event nb_event] [-l lvllog] [--log lvllog] [-h]  OutputFile RootInputFile_withParObj [RootInputFiles....]"<<std::endl; std::exit(1); break;
      case 'c': std::cout<<"Nb CPU "<< optarg<<std::endl; Nb_CPU=std::atoi(optarg); break;
      case 'n': std::cout<<"fraction of event "<<optarg<<std::endl; Nb_fraction=std::atoi(optarg); break;
      case 'e': std::cout<<"Nb Event "<<optarg<<std::endl; Nb_event=std::atoi(optarg); break;
      case 's': std::cout<<"Start Event "<<optarg<<std::endl; Start=std::atoi(optarg); break;
      case 'g': std::cout<<"Geometry field :"<<optarg<<std::endl; nameGeo= std::string(optarg); break;
      case 'l': std::cout<<"Log level:"<<optarg<<std::endl; Log_lvl = std::atoi(optarg); break;
      case '?': std::cerr <<"usage: "<<argv[0]<<" [-g Geofile] [--geo Geofile] [-c nb_cpu] [--cpu nb_cpu] [-n fraction] [--num fraction] [-s start_ev] [--start start_ev] [-e nb_event] [--event nb_event] [-l lvllog] [--log lvllog] [-h]  OutputFile RootInputFile "<<std::endl; std::exit(1);
      }
	
  std::string name_in,name_out;

  if(optind == argc )
    {
		
      std::cerr <<"usage: "<<argv[0]<<" [-g Geofile] [--geo Geofile] [-c nb_cpu] [--cpu nb_cpu] [-n fraction] [--num fraction] [-s start_ev] [--start start_ev] [-e nb_event] [--event nb_event] [-l lvllog] [--log lvllog] [-h]  OutputFile RootInputFile "<<std::endl;
      std::cerr <<" input and output Rootfile are missing !"<<optind<<" "<<argc<<std::endl;
      std::exit(1);
    }
  else
    {
      name_out = argv[optind];
      //if(EventGenConf==false)
      name_in = argv[optind+1];
    }

  
  try
    {
      auto sink1 = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
      sink1->set_pattern("[%^%l%$] %v");
      auto ConsoleLogger = std::make_shared<spdlog::logger>("Console", sink1);
      spdlog::register_logger(ConsoleLogger);
      //auto console = spdlog::stdout_color_mt("Console");
      switch (Log_lvl)
	{
	case -1:
	  ConsoleLogger->set_level(spdlog::level::off); break;
	case 0:
	  ConsoleLogger->set_level(spdlog::level::warn); break;
	case 1:
	  ConsoleLogger->set_level(spdlog::level::info); break;
	case 2:
	  ConsoleLogger->set_level(spdlog::level::debug); break;
	default: 
	  ConsoleLogger->set_level(spdlog::level::info); break;
	}
      
      std::string file_base_name = name_out.substr(name_out.find_last_of('/')+1);
      ConsoleLogger->info("name in : {} / name out : {}", name_in, name_out);

      //  std::cout<<"name in :"<<name_in<<" / name out :"<<name_out<<std::endl;
  
      //TFile *sourcefile = new TFile (name_out_source.c_str(),"RECREATE");

    
      //sourcefile->Close();
      //delete sourcefile;
  

      /*******************************************************************/

      TStopwatch timer;
      timer.Start(kFALSE);

      // HypHiMC_output_TREE ***********************
      Long64_t total_nentries=0;
      //Long64_t total_nentries_over_files=0;
      int Start_event[3]={0,0,0};
      int Stop_event[3]={0,0,0};



      TGeoManager::Import(nameGeo.c_str()); // FOR Kalman_DAF

      TDatabasePDG::Instance()->AddParticle("deuteron","deuteron",1.875613,kTRUE,0.,1.*3.,"Ions",10000);
      TDatabasePDG::Instance()->AddParticle("triton","triton",2.80892/*2.80925*/,kTRUE,0.,1.*3.,"Ions",10001);
      TDatabasePDG::Instance()->AddParticle("alpha","alpha",3.72738/*3.727417*/,kTRUE,0.,2.*3.,"Ions",10002);
      TDatabasePDG::Instance()->AddParticle("He3","He3",2.80839/*2.80923*/,kTRUE,0.,2.*3.,"Ions",10003);
      TDatabasePDG::Instance()->AddParticle("Li6","Li6",5.60152/*5.6015194*/,kTRUE,0.,3.*3.,"Ions",10004);
      TDatabasePDG::Instance()->AddParticle("H3L","H3L",2.99114,kFALSE,0.,1.*3.,"Ions",20001);
      TDatabasePDG::Instance()->AddParticle("H4L","H4L",3.9225,kFALSE,0.,1.*3.,"Ions",20002);
      TDatabasePDG::Instance()->AddParticle("He5L","He5L",4.8399,kFALSE,0.,1.*3.,"Ions",20003);
  

      TFile* offile = new TFile(name_out.c_str(),"RECREATE");

      TTree* Anatree = 0;
      MCAnaEventG4Sol *Anaevent = 0;
      Anatree = new TTree("T","MCAnaEventG4SolTree");
      //Anatree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
      Anatree->SetAutoSave(100000000); // autosave when 100 Mbyte written
      //Anatree->SetAutoSave(50000000); // autosave when 50 Mbyte written
      //Anatree->SetMaxTreeSize(50000000); //1 GB
      Anatree->SetAutoFlush();
      Anaevent = new MCAnaEventG4Sol();
      //Anatree->Bronch("bAna_event","MCAnaEvent",&Anaevent, 32000,1);
      MCAnaEventG4Sol::Class()->IgnoreTObjectStreamer();
      TBranch* branch = Anatree->Branch("MCAnaEventG4Sol",&Anaevent, 32000,2);
      branch->SetAutoDelete(kFALSE);
      Anatree->BranchRef();

      long long int nb =0;

      TFile* input_file = new TFile(name_in.c_str());

      TTree* InTree = 0;
      InTree = dynamic_cast<TTree*> (input_file->Get("G4Tree"));
      assert(InTree!=0);

      DataSim InputPar{nullptr,nullptr}; 
      InputPar.nameDet = (std::vector<std::string>*)(input_file->Get("nameDet"));
      InputPar.simParameters = (std::map<std::string,double>*)(input_file->Get("simParameters"));
  
      if(InputPar.nameDet == nullptr)
	{
	  ConsoleLogger->error("E> Load nameDet not possible ! {}", fmt::ptr(input_file->Get("nameDet")));
	  return -1;
	}
      if(InputPar.simParameters == nullptr)
	{
	  ConsoleLogger->error("E> Load simParameters not possible ! {}", fmt::ptr(input_file->Get("simParameters")));
	  return -1;
	}
      TG4Sol_Event* fEvent = 0;
      InTree->SetBranchAddress("TG4Sol_Event", &fEvent);

#ifdef ROOT6
      TTreeReader reader(InTree);

      TTreeReaderValue<TG4Sol_Event> ReaderEvent(reader,"TG4Sol_Event");
      std::vector<TTreeReaderArray<TG4Sol_Hit>*> AllHits;
      for(auto name : *InputPar.nameDet)
	AllHits.emplace_back(new TTreeReaderArray<TG4Sol_Hit>(reader,name.c_str()));
#else
      std::vector<TClonesArray*> AllHits;
      for(auto name : *InputPar.nameDet)
	{
	  AllHits.emplace_back( new TClonesArray("TG4Sol_Hit",20));
	  InTree->SetBranchAddress(name.c_str(),&AllHits.back());
	  AllHits.back()->SetName(name.c_str());
	}
#endif 

#ifdef TREEPERF
      TFile* f_in = new TFile(argv[optind+1]);
      TTree* Chain_InTree = (TTree*) f_in->Get("HypHiMC_output_TREE");
#endif

      Long64_t total_nentries_over_files = InTree->GetEntries();
      int max_loop = 1;
  
      Start_event[0]=Start;
      Stop_event[0]=Nb_event !=-1 ? Nb_event+Start : total_nentries_over_files;
    
      if(Stop_event[0]>total_nentries_over_files)
	Stop_event[0]=total_nentries_over_files;
      if(Start_event[0]>total_nentries_over_files)
	Start_event[0]=total_nentries_over_files;

      std::list<std::string> Option;
      std::list<std::string> ToDo;
      TString cpu_string("CPU");
      cpu_string+=Nb_CPU;
      Option.push_back(cpu_string.Data());
      TString fra_string("Fraction");
      fra_string+=Nb_fraction;
      Option.push_back(fra_string.Data());
      TString nevent_string("Event");
      nevent_string+=Stop_event[0]-Start_event[0];
      Option.push_back(nevent_string.Data());
      //Option.push_back("NoMaterial");
      Option.push_back("Debug_DAF");
      
      
      Ana_Hist ListHisto(true/*DAF_Debug*/,false/*Oldvertex*/,false/*DCproject*/,true/*Finding*/,true/*Hough*/,true);
      FullRecoTask ReconstructionTask(ToDo,Option,1.,InputPar);

      ReconstructionTask.AttachHisto(&ListHisto);
      ConsoleLogger->info("Init done");
      //-------------------------------------------------------
      // EVENT LOOP nentries
      TTimeStamp stamp;

      ConsoleLogger->info("Start = {} Stop= {}", Start_event[0], Stop_event[0]);
      ConsoleLogger->info("Start = {} Stop= {}", Start_event[1], Stop_event[1]);
      ConsoleLogger->info("Start = {} Stop= {}", Start_event[2], Stop_event[2]);
      ConsoleLogger->info("---------------------------------------------");

      Long64_t for_display = 0;
      Long64_t timing = 0;
      //reader.SetEntriesRange(Start_event[0],Stop_event[0]);

      for(Long64_t i = Start_event[0];i<Stop_event[0];++i)
	{
	  InTree->GetEntry(i);

	  //while(reader.Next())
	  //{  
	  Long64_t iEvent = i;//reader.GetCurrentEntry(); 
      
	  ++total_nentries;
	  Anaevent->Clear();
      
	  if(iEvent%10000==0)
	    ConsoleLogger->info("Processing Event# {} / {} | {} %", iEvent, Stop_event[0], (double)iEvent/(double)(Stop_event[0]-Start_event[0])*100);
	  
	  if((int)((double)iEvent/(double)(Stop_event[0]-Start_event[0])*4)==timing)
	    {
	      ConsoleLogger->info("Progress : {} %", (int)((double)iEvent/(double)(Stop_event[0]-Start_event[0])*100.));
	      ++timing; 
	    }
	  //auto event = ReaderEvent.Get();
      
	  int toStop = ReconstructionTask.EventLoop(*fEvent, AllHits, Anaevent);
      
	  nb += Anatree->Fill();
	  ++for_display;
	  Anaevent->Clear();
  
	}
  

      offile=Anatree->GetCurrentFile();
      offile->cd();
      Anatree->Write(); 
      ConsoleLogger->info("Tree Filled");
      //ListHisto.Write(offile_hist);
      ListHisto.Write(offile);

      ConsoleLogger->info("Histo Written");

      TObjArray* additionalPDG = new TObjArray;
      additionalPDG->SetName("additionalPDG");
      int id_pdg = 10004+1;
      bool finish = false;
      while( finish==false )
	{
	  TParticlePDG* PDG_particle = TDatabasePDG::Instance()->GetParticle(id_pdg);
	  if(PDG_particle!=nullptr)
	    {
	      additionalPDG->Add(PDG_particle);
	      ++id_pdg;
	    }
	  else
	    finish = true;
	}
      id_pdg = 20003+1;
      while( finish==false )
	{
	  TParticlePDG* PDG_particle = TDatabasePDG::Instance()->GetParticle(id_pdg);
	  if(PDG_particle!=nullptr)
	    {
	      additionalPDG->Add(PDG_particle);
	      ++id_pdg;
	    }
	  else
	    finish = true;
	}
      offile->WriteObjectAny(additionalPDG, additionalPDG->Class(), "additionalPDG");
      //additionalPDG->Write();

      ConsoleLogger->info("Additional PDG data written");

      offile->Close();
      //offile_hist->Close();
      ConsoleLogger->info("Histo File Closed");

#ifdef TREEPERF
      ps->SaveAs("perfstat_2.root");
#endif
      

      delete offile; offile=NULL;
      ConsoleLogger->info("Offile deleted");
      // delete offile_hist;offile_hist=NULL;
      //std::cout<<"Offile Histo deleted"<<std::endl;
  
      delete Anaevent; Anaevent=NULL;
      ConsoleLogger->info("Anaevent deleted");

#ifdef TREEPERF
      delete ps; ps=0;
#endif
      /***************************************************************/
      timer.Stop();
      //delete AnaOutTree;
  

      Float_t mbytes = 0.000001*nb;
      Double_t rtime = timer.RealTime();
      Double_t ctime = timer.CpuTime();
      ConsoleLogger->info("");
      ConsoleLogger->info("{} events and {} bytes processed.", total_nentries ,nb);
      ConsoleLogger->info("RealTime= {} seconds, CpuTime= {} seconds", rtime, ctime);
      ConsoleLogger->info("You read {} Mbytes/Realtime seconds", mbytes/rtime);
      ConsoleLogger->info("You read {} Mbytes/Cputime seconds", mbytes/ctime);

      // Saving histograms
  
      //AnaOutTree->Write(AnaOutTree->name_outfile);

      //delete Anaevent;
      //Anatree->Reset();
      //delete Anatree;

      //delete app;

      return 0;
    }
  catch (const spdlog::spdlog_ex& ex)
    {
      std::cout << "Log initialization failed: " << ex.what() << std::endl;
      std::exit(-1);
    }
}

