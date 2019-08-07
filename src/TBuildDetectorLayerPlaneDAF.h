#ifndef TDATABUILDLPDAF
#define TDATABUILDLPDAF

#include "TDataBuilder.h"

#include "THyphiAttributes.h"
#include "TRandom3.h"

#include "FullRecoEvent.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "EventG4Sol/TG4Sol_Event.hh"
#include "EventG4Sol/TG4Sol_Hit.hh"

//#include "MathematicalTools.hh"
#include "Debug.hh"

#include "TProfile.h"
#include "TFile.h"

#include "TGeoManager.h"
#include "TGeoElement.h"

#include <sstream>

#include "spdlog/spdlog.h"

//#include "TClustering.h"
#ifdef ROOT6
#include "TTreeReaderArray.h"
#endif

class PDG_fromName {

  static const std::string ElName2[];
  
  TGeoElementTable *tableRN; 
  std::unordered_map<std::string,double> cache_NucleiPID;

  const double u = 0.931494061;
  const double m_eminus = 5.10998909999999971e-04; // GeV

  int lastPDGIon = 10004; 
public:
  PDG_fromName():tableRN(nullptr)
  {  }
  ~PDG_fromName()
  {
    std::stringstream ss;
    ss<<"PDG_fromName cache \n";
    for(auto it : cache_NucleiPID)
      ss<<"["<<it.first<<", "<<it.second<<"] ";
    spdlog::get("Console")->info(ss.str());
  }
  
  int operator() (const std::string& name)
  {
    //std::cout<<"!> PID from name:"<<name<<"\n";
    
    if(name=="proton")
      return 2212;
    if(name=="He3")
      return 10003;
    if(name=="pi-")
      return -211;
    if(name=="pi+")
      return 211;
    if(name=="neutron")
      return 2112;
    if(name=="kaon0")
      return 311;
    if(name=="kaon+")
      return 321;
    if(name=="kaon-")
      return -321;
    if(name=="H3L")
      return 20001;
    if(name=="pi0")
      return 111;
    if(name=="triton")
      return 10001;
    if(name=="deuteron")
      return 10000;
    if(name=="alpha")
      return 10002;
    if(name=="He3")
      return 10003;
    if(name=="Li6")
      return 10004;

    
    auto it_name = cache_NucleiPID.find(name);
    if(it_name!=cache_NucleiPID.end())
      return it_name->second;

    //std::cout<<"ADD NEW ELEM TO PDG DATABASE: "<<name<<"\n";
    
    std::size_t posCharge = name.find_first_of("0123456789");

    if(posCharge == std::string::npos)
      return 0;

    //std::cout<<"found Charge at"<<posCharge;
    int AtomMass = std::stoi(name.substr(posCharge,std::string::npos));
    std::string nameElement = name.substr(0,posCharge);

    //std::cout<<" Unpack into:"<<nameElement<<" "<<AtomMass<<"\n";

    int id_Elem = -1;
    for(size_t i=0; i<111; ++i)
      if(nameElement == ElName2[i])
	{
	  id_Elem=i;
	  break;
	}
    if(id_Elem<=0)
      return 0;

    //std::cout<<"Element found :"<<id_Elem<<" "<<ElName2[id_Elem]<<"\n";
    
    if(tableRN==nullptr)
      tableRN = gGeoManager->GetElementTable();    
    
    auto* TempElement = tableRN->GetElementRN(AtomMass,id_Elem);
    double Dmass = 0.;
    if(TempElement==nullptr)
      {
	spdlog::get("some_logger")->info("E> no element ! {} {} {}", AtomMass, id_Elem, fmt::ptr(TempElement));
	if(AtomMass==23 && id_Elem==14)
	  Dmass = 23.073*1e-3; // MeV -> GeV
	else
	  return 0;
      }
    else
      Dmass = TempElement->MassEx()*1e-3; // MeV -> GeV
    
    double Mass = Dmass+AtomMass*u - id_Elem*m_eminus;
    
    TDatabasePDG::Instance()->AddParticle(name.c_str(),name.c_str(),Mass,kFALSE,0.,id_Elem*3.,"Ions",++lastPDGIon);

    cache_NucleiPID.insert(std::make_pair(name,lastPDGIon));

    return lastPDGIon;
    
    //return 0;
  }

};

constexpr bool IsPlanar(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::InSi0 : ;
  case G4Sol::InSi1 : ;
  case G4Sol::InSi2 : ;
  case G4Sol::InSi3 : ;
  case G4Sol::TR1 : ;
  case G4Sol::TR2 : ;
  case G4Sol::PSFE : ;
  case G4Sol::PSBE : ;
  case G4Sol::TrFwd0 : ;
  case G4Sol::TrFwd1 : ;
  case G4Sol::TrFwd2 : ;
  case G4Sol::RPC_l : ;
  case G4Sol::RPC_h : ;
  case G4Sol::FMF2Stop0 : ;
  case G4Sol::FMF2Stop1 : ;
  case G4Sol::FMF2Stop2 :
    return true;
  default:
    return false;
  };  
}; 



class TBuildDetectorLayerPlaneDAF : public TDataBuilder
{

 public:
  
  const THyphiAttributes& att;

  TBuildDetectorLayerPlaneDAF(const THyphiAttributes& att);
  ~TBuildDetectorLayerPlaneDAF();

#ifdef ROOT6
  int operator() (const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#else
  int operator() (const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#endif
  
  private :
  //int Exec(THyphiEvent_Phys_new *event,Ana_Event* OutTree);
#ifdef ROOT6
  int Exec(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#else  
  int Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#endif
  
  int SoftExit(int);

private:
  
  //std::vector<std::vector<GFAbsRecoHit*> > ListAllHits;
  std::unordered_map<int,int> orderDetectors;
  PDG_fromName pid_fromName;
};

#endif
