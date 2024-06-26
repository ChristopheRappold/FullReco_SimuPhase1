#ifndef TDATABUILDLPDAF_MT
#define TDATABUILDLPDAF_MT

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "EventG4Sol/TG4Sol_Event.hh"
#include "EventG4Sol/TG4Sol_Hit.hh"
#include "FullRecoEvent.hh"
#include "TDataBuilder.h"
#include "THyphiAttributes.h"
#include "TRandom3.h"

//#include "MathematicalTools.hh"
#include "Debug.hh"
#include "TFile.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TProfile.h"
#include "spdlog/spdlog.h"

#include <sstream>

//#include "TClustering.h"
#ifdef ROOT6
#include "TTreeReaderArray.h"
#endif

// class PDG_fromName
// {

//   static const std::string ElName2[];

//   TGeoElementTable* tableRN;
//   std::unordered_map<std::string, double> cache_NucleiPID;

//   const double u        = 0.931494061;
//   const double m_eminus = 5.10998909999999971e-04; // GeV

//   int lastPDGIon = 10004;

// public:
//   PDG_fromName() : tableRN(nullptr) {}
//   ~PDG_fromName()
//   {
//     std::stringstream ss;
//     ss << "PDG_fromName cache \n";
//     for(const auto& it : cache_NucleiPID)
//       ss << "[" << it.first << ", " << it.second << "] ";
//     spdlog::get("Console")->info(ss.str());
//   }

//   int operator()(const std::string& name)
//   {
//     // std::cout<<"!> PID from name:"<<name<<"\n";

//     if(name == "proton")
//       return 2212;
//     if(name == "He3")
//       return 10003;
//     if(name == "pi-")
//       return -211;
//     if(name == "pi+")
//       return 211;
//     if(name == "neutron")
//       return 2112;
//     if(name == "kaon0")
//       return 311;
//     if(name == "kaon+")
//       return 321;
//     if(name == "kaon-")
//       return -321;
//     if(name == "H3L")
//       return 20001;
//     if(name == "pi0")
//       return 111;
//     if(name == "triton")
//       return 10001;
//     if(name == "deuteron")
//       return 10000;
//     if(name == "alpha")
//       return 10002;
//     if(name == "He3")
//       return 10003;
//     if(name == "Li6")
//       return 10004;

//     auto it_name = cache_NucleiPID.find(name);
//     if(it_name != cache_NucleiPID.end())
//       return it_name->second;

//     // std::cout<<"ADD NEW ELEM TO PDG DATABASE: "<<name<<"\n";

//     std::size_t posCharge = name.find_first_of("0123456789");

//     if(posCharge == std::string::npos)
//       return 0;

//     // std::cout<<"found Charge at"<<posCharge;
//     int AtomMass            = std::stoi(name.substr(posCharge, std::string::npos));
//     std::string nameElement = name.substr(0, posCharge);

//     // std::cout<<" Unpack into:"<<nameElement<<" "<<AtomMass<<"\n";

//     int id_Elem = -1;
//     for(size_t i = 0; i < 111; ++i)
//       if(nameElement == ElName2[i])
//         {
//           id_Elem = i;
//           break;
//         }
//     if(id_Elem <= 0)
//       return 0;

//     // std::cout<<"Element found :"<<id_Elem<<" "<<ElName2[id_Elem]<<"\n";

//     if(tableRN == nullptr)
//       tableRN = gGeoManager->GetElementTable();

//     auto* TempElement = tableRN->GetElementRN(AtomMass, id_Elem);
//     double Dmass      = 0.;
//     if(TempElement == nullptr)
//       {
//         spdlog::get("some_logger")->info("E> no element ! {} {} {}", AtomMass, id_Elem, fmt::ptr(TempElement));
//         if(AtomMass == 23 && id_Elem == 14)
//           Dmass = 23.073 * 1e-3; // MeV -> GeV
//         else
//           return 0;
//       }
//     else
//       Dmass = TempElement->MassEx() * 1e-3; // MeV -> GeV

//     double Mass = Dmass + AtomMass * u - id_Elem * m_eminus;

//     TDatabasePDG::Instance()->AddParticle(name.c_str(), name.c_str(), Mass, kFALSE, 0., id_Elem * 3., "Ions",
//                                           ++lastPDGIon);

//     cache_NucleiPID.insert(std::make_pair(name, lastPDGIon));

//     return lastPDGIon;

//     // return 0;
//   }
// };

constexpr bool IsPlanar(G4Sol::SolDet idDet)
{
  switch(idDet)
    {
    case G4Sol::InSi0:;
    case G4Sol::InSi1:;
    case G4Sol::InSi2:;
    case G4Sol::InSi3:;
    case G4Sol::TR1:;
    case G4Sol::TR2:;
    case G4Sol::PSFE:;
    case G4Sol::PSBE:;
    case G4Sol::TrFwd0:;
    case G4Sol::TrFwd1:;
    case G4Sol::TrFwd2:;
    case G4Sol::RPC_l:;
    case G4Sol::RPC_h:;
    case G4Sol::FMF2Stop0:;
    case G4Sol::FMF2Stop1:;
    case G4Sol::FMF2Stop2:
      return true;
    default:
      return false;
    };
};

constexpr bool IsPSCE(G4Sol::SolDet idDet)
{
  switch(idDet)
    {
    case G4Sol::PSCE:
      return true;
    default:
      return false;
    };
};
constexpr bool IsFiberU_Vetoed(G4Sol::SolDet idDet)
{
  switch(idDet)
    {
    case G4Sol::FiberD1_x:;
    case G4Sol::FiberD1_u:;
    case G4Sol::FiberD1_v:;
    case G4Sol::FiberD2_x:;
    case G4Sol::FiberD2_u:;
    case G4Sol::FiberD2_v:;
    case G4Sol::FiberD3_x:;
    case G4Sol::FiberD3_u:;
    case G4Sol::FiberD3_v:;
      return true;
    default:
      return false;
    };
};

constexpr bool IsFiberU(G4Sol::SolDet idDet)
{
  switch(idDet)
    {
    case G4Sol::FiberD1_x:;
    case G4Sol::FiberD1_u:;
    case G4Sol::FiberD1_v:;
    case G4Sol::FiberD2_x:;
    case G4Sol::FiberD2_u:;
    case G4Sol::FiberD2_v:;
    case G4Sol::FiberD3_x:;
    case G4Sol::FiberD3_u:;
    case G4Sol::FiberD3_v:;
    case G4Sol::FiberD4_v:;
    case G4Sol::FiberD4_u:;
    case G4Sol::FiberD4_x:;
    case G4Sol::FiberD5_x:;
    case G4Sol::FiberD5_u:;
    case G4Sol::FiberD5_v:;
    case G4Sol::MiniFiberD1_x:;
    case G4Sol::MiniFiberD1_v:;
    case G4Sol::MiniFiberD1_u:;
    case G4Sol::MiniFiberD2_x:;
    case G4Sol::MiniFiberD2_u:;
    case G4Sol::MiniFiberD2_v:;
      return true;
    default:
      return false;
    };
};

constexpr bool IsWire(G4Sol::SolDet idDet)
{
  switch(idDet)
    {
    case G4Sol::MG01:;
    case G4Sol::MG02:;
    case G4Sol::MG03:;
    case G4Sol::MG04:;
    case G4Sol::MG05:;
    case G4Sol::MG06:;
    case G4Sol::MG07:;
    case G4Sol::MG08:;
    case G4Sol::MG09:;
    case G4Sol::MG10:;
    case G4Sol::MG11:;
    case G4Sol::MG12:;
    case G4Sol::MG13:;
    case G4Sol::MG14:;
    case G4Sol::MG15:;
    case G4Sol::MG16:;
    case G4Sol::MG17:;
      return true;
    default:
      return false;
    };
};

double CloseDist(const TVector3& Xin, const TVector3& Xout, const TVector3& Pin, const TVector3& Pout);

class TBuildDetectorLayerPlaneDAF_MT final : public TDataBuilder
{

public:
  const THyphiAttributes& att;

  explicit TBuildDetectorLayerPlaneDAF_MT(const THyphiAttributes& att);
  ~TBuildDetectorLayerPlaneDAF_MT() final;

#ifdef ROOT6
  ReturnRes::InfoM operator()(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                              FullRecoEvent& RecoEvent);
#else
  ReturnRes::InfoM operator()(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits,
                              FullRecoEvent& RecoEvent);
#endif

  double CloseDist(const TVector3& Xin, const TVector3& Xout, const TVector3& Pin, const TVector3& Pout);

private:
  // int Exec(THyphiEvent_Phys_new *event,Ana_Event* OutTree);
#ifdef ROOT6
  int Exec(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, FullRecoEvent& RecoEvent);
#else
  int Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent);
#endif

  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

private:
  // std::vector<std::vector<GFAbsRecoHit*> > ListAllHits;
  std::unordered_map<int, int> orderDetectors;
  std::unordered_map<int, std::string> orderDetName;
  PDG_fromName pid_fromName;
  struct LocalHists
  {
    TH1I* h_stats;
  };
  LocalHists LocalHisto;
};

#endif
