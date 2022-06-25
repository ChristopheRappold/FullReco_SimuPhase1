#ifndef TDATABUILDRESTART
#define TDATABUILDRESTART

#include "TDataBuilder.h"

#include "THyphiAttributes.h"
#include "TRandom3.h"

#include "FullRecoEvent.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "Debug.hh"

#include "TProfile.h"
#include "TFile.h"

#include "TGeoManager.h"
#include "TGeoElement.h"

#include <sstream>
#include <type_traits>

#include "spdlog/spdlog.h"

//#include "TClustering.h"
#ifdef ROOT6
#include "TTreeReaderArray.h"
#endif

namespace restart {

constexpr bool IsFiberU_Vetoed(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::FiberD1_x : ;
  case G4Sol::FiberD1_u : ;
  case G4Sol::FiberD1_v : ;
  case G4Sol::FiberD2_x : ;
  case G4Sol::FiberD2_u : ;
  case G4Sol::FiberD2_v : ;
  case G4Sol::FiberD3_x : ;
  case G4Sol::FiberD3_u : ;
  case G4Sol::FiberD3_v : ;
    return true;
  default:
    return false;
  };
};

constexpr bool IsFiberU(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::FiberD1_x : ;
  case G4Sol::FiberD1_u : ;
  case G4Sol::FiberD1_v : ;
  case G4Sol::FiberD2_x : ;
  case G4Sol::FiberD2_u : ;
  case G4Sol::FiberD2_v : ;
  case G4Sol::FiberD3_x : ;
  case G4Sol::FiberD3_u : ;
  case G4Sol::FiberD3_v : ;
  case G4Sol::FiberD4_x : ;
  case G4Sol::FiberD4_u : ;
  case G4Sol::FiberD4_v : ;
  case G4Sol::FiberD5_x : ;
  case G4Sol::FiberD5_u : ;
  case G4Sol::FiberD5_v : ;
  // case G4Sol::MiniFiberD1_x1 : ;
  // case G4Sol::MiniFiberD1_u1 : ;
  // case G4Sol::MiniFiberD1_v1 : ;
  // case G4Sol::MiniFiberD1_x2 : ;
  // case G4Sol::MiniFiberD1_u2 : ;
  // case G4Sol::MiniFiberD1_v2 : ;
    return true;
  default:
    return false;
  };
};


constexpr bool IsFiberM(G4Sol::SolDet idDet) {
  switch (idDet) {
    case G4Sol::MiniFiberD1_x1 : ;
    case G4Sol::MiniFiberD1_u1 : ;
    case G4Sol::MiniFiberD1_v1 : ;
    case G4Sol::MiniFiberD1_x2 : ;
    case G4Sol::MiniFiberD1_u2 : ;
    case G4Sol::MiniFiberD1_v2 : ;
                                 return true;
    default:
                                 return false;
  };
};

double CloseDist( const TVector3 & Xin, const TVector3 & Xout,
		  const TVector3 & Pin, const TVector3 & Pout );


// C++20 only
// bool has_fMC_Particle(const auto &obj)
// {
//   if constexpr (requires {obj.fMC_Particle;})
//     {
//       return true;
//     }
//   else
//       return false;
// };


template <typename T, typename = int>
struct HasMC_Particle : std::false_type { };

template <typename T>
struct HasMC_Particle <T, decltype((void) T::fMC_Particle, 0)> : std::true_type { };

}
template<typename Out>
class TBuildRestarter final : public TDataBuilder
{

public:

  const THyphiAttributes& att;

  explicit TBuildRestarter(const THyphiAttributes& att);
  ~TBuildRestarter() final;

#ifdef ROOT6
  ReturnRes::InfoM operator() (const Out& RestartEvent, FullRecoEvent& RecoEvent, Out* OutTree);
#else
  ReturnRes::InfoM operator() (Out* event, FullRecoEvent& RecoEvent, Out* OutTree);
#endif

private :
  //int Exec(THyphiEvent_Phys_new *event,Ana_Event* OutTree);
#ifdef ROOT6
  int Exec(const Out& RestartEvent, FullRecoEvent& RecoEvent, Out* OutTree);
#else
  int Exec(Out* event, FullRecoEvent& RecoEvent, Out* OutTree);
#endif

  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

private:

  //std::vector<std::vector<GFAbsRecoHit*> > ListAllHits;
  std::unordered_map<int,int> orderDetectors;
  std::unordered_map<int,std::string> orderDetName;
  PDG_fromName pid_fromName;
  struct LocalHists
  {
    TH1I* h_stats;
  };
  LocalHists LocalHisto;
};

#endif
