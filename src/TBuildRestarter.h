#ifndef TDATABUILDRESTART
#define TDATABUILDRESTART

#include "TDataBuilder.h"

#include "THyphiAttributes.h"
#include "TRandom3.h"

#include "FullRecoEvent.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "Ana_Event/Ana_WasaEvent.hh"

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

}

class TBuildRestarter final : public TDataBuilder
{

public:

  const THyphiAttributes& att;

  explicit TBuildRestarter(const THyphiAttributes& att);
  ~TBuildRestarter() final;

#ifdef ROOT6
  ReturnRes::InfoM operator() (const MCAnaEventG4Sol& RestartEvent, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#else
  ReturnRes::InfoM operator() (MCAnaEventG4Sol* event, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#endif


#ifdef ROOT6
  ReturnRes::InfoM operator() (const Ana_WasaEvent& RestartEvent, FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree) { }
#else
  ReturnRes::InfoM operator() (Ana_WasaEvent* event, FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree) { }
#endif


private :
  //int Exec(THyphiEvent_Phys_new *event,Ana_Event* OutTree);
#ifdef ROOT6
  int Exec(const MCAnaEventG4Sol& RestartEvent, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#else
  int Exec(MCAnaEventG4Sol* event, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
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
