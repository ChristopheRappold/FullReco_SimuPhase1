#ifndef TDATABUILDLPDAF
#define TDATABUILDLPDAF

#include "TDataBuilder.h"

#include "THyphiAttributes.h"
#include "TRandom3.h"

#include "FullRecoEvent.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "EventG4Sol/TG4Sol_Event.hh"
#include "EventG4Sol/TG4Sol_Hit.hh"

#include "Ana_Event/Ana_WasaEvent.hh"
#include "HitAna/FiberAnalyzer.hh"
#include "HitAna/FiberHitAna.hh"
#include "HitAna/FiberHitXUV.hh"
#include "HitAna/FiberTrackAna.hh"

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


constexpr bool IsPlanar(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::InSi0 : ;
  case G4Sol::InSi1 : ;
  case G4Sol::InSi2 : ;
  case G4Sol::InSi3 : ;
  case G4Sol::TR1 : ;
  case G4Sol::TR2 : ;
  case G4Sol::Si1x: ;
  case G4Sol::Si1y: ;
  case G4Sol::Si2x: ;
  case G4Sol::Si2y: ;
  case G4Sol::Si1x_SD: ;
  case G4Sol::Si1y_SD: ;
  case G4Sol::Si2x_SD: ;
  case G4Sol::Si2y_SD: ;
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

constexpr bool IsPSCE(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::PSCE :
    return true;
  default:
    return false;
  };
};
constexpr bool IsPSBE(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::PSBE :
    return true;
  default:
    return false;
  };
};

constexpr bool IsSilicon(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::Si1x : ;
  case G4Sol::Si1y : ;
  case G4Sol::Si2x : ;
  case G4Sol::Si2y : ;
    return true;
  default:
    return false;
  };
};

constexpr bool IsSilicon_SD(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::Si1x_SD : ;
  case G4Sol::Si1y_SD : ;
  case G4Sol::Si2x_SD : ;
  case G4Sol::Si2y_SD : ;
    return true;
  default:
    return false;
  };
};

constexpr bool IsSilicon_SD_pad(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::Si1x_SD_pad : ;
  case G4Sol::Si1y_SD_pad : ;
  case G4Sol::Si2x_SD_pad : ;
  case G4Sol::Si2y_SD_pad : ;
    return true;
  default:
    return false;
  };
};

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

constexpr bool IsWire(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::MG01 : ;
  case G4Sol::MG02 : ;
  case G4Sol::MG03 : ;
  case G4Sol::MG04 : ;
  case G4Sol::MG05 : ;
  case G4Sol::MG06 : ;
  case G4Sol::MG07 : ;
  case G4Sol::MG08 : ;
  case G4Sol::MG09 : ;
  case G4Sol::MG10 : ;
  case G4Sol::MG11 : ;
  case G4Sol::MG12 : ;
  case G4Sol::MG13 : ;
  case G4Sol::MG14 : ;
  case G4Sol::MG15 : ;
  case G4Sol::MG16 : ;
  case G4Sol::MG17 : ;
    return true;
  default:
    return false;
  };
};

double CloseDist( const TVector3 & Xin, const TVector3 & Xout,
		  const TVector3 & Pin, const TVector3 & Pout );




class TBuildDetectorLayerPlaneDAF final : public TDataBuilder
{

public:

  const THyphiAttributes& att;

  explicit TBuildDetectorLayerPlaneDAF(const THyphiAttributes& att);
  ~TBuildDetectorLayerPlaneDAF() final;

#ifdef ROOT6
  ReturnRes::InfoM operator() (const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#else
  ReturnRes::InfoM operator() (const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#endif

ReturnRes::InfoM operator()(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree);

private :
  //int Exec(THyphiEvent_Phys_new *event,Ana_Event* OutTree);
#ifdef ROOT6
  int Exec(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#else
  int Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#endif
  
  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

private:
  
  //std::vector<std::vector<GFAbsRecoHit*> > ListAllHits;
  std::unordered_map<int,int> orderDetectors;
  std::unordered_map<int,std::string> orderDetName;
  PDG_fromName pid_fromName;
  int offsetGeoNameID_MDC = 0;
  int offsetGeoNameID_PSCE = 0;
  struct LocalHists
  {
    TH1I* h_stats;
  };
  
  LocalHists LocalHisto;

public:

  std::unique_ptr<ParaManager> par;

};

#endif
